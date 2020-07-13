/* The multiple-reactivation model with the first recrudescence time
 * as a latent variable.
 * The time of treatment initiation is added as a covariate
 * for the recrudescence rate lambda and the growth rate r (a.k.a. g).
 * 
 * Abbreviations:
 * Cndl: conditional
 * Fst: first
 * FPT: first passage time
 * Pts: points
 */


/* In the functions block, we can define some custom funtions and distribution.
 * In our case, we define the rebound-time distribution f(t;...) (called fpt below)
 * and some other auxiliary functions
 */
functions {
    /* A "potentially censored" normal distribution.
     * Typically used for virus load measurements 
     * that can fall below a dectection limit.
     * The type of censoring is determined by the control
     * parameter cc.
     */
    real censored_normal_lpdf(real x, real xhat, real sigma, int cc) {
        real lp;
        if ( cc == 0 ) {
            lp = normal_lpdf(x | xhat, sigma);
        } else if ( cc == 1 ) {
            lp = normal_lcdf(x | xhat, sigma);
        } else if ( cc == 2 ) {
            lp = normal_lccdf(x | xhat, sigma);
        } else if ( cc == 3 ) {
            lp = 0;
        } else {
            reject("invalid censor code");
        }
        return lp;
    }
    /* Compute the log of the logistic growth curve at time points ts,
     * given parameters r, logK, logTheta, log1mTheta and tm.
     * Here theta = ell / K, and V(tm) = ell, and 
     * log1mTheta = log(1 - theta).
     * We use the log_sum_exp function for numeric stability.
     */
    vector log_logistic_model(vector ts, real r, real logK, 
            real logTheta, real log1mTheta, real tm) {
        vector[num_elements(ts)] xs;
        for ( j in 1:num_elements(ts) ) {
            xs[j] =  logTheta + logK - log_sum_exp(-r*(ts[j]-tm) + log1mTheta, logTheta);
        }
        return xs;
    }
    /* A simple function to create a sequence of n equally 
     * spaced numbers between xmin and xmax (inclusive).
     * This is used in the generated quantities block for simulations.
     */
    vector seq(real xmin, real xmax, int n) { // n should be > 2
        real mesh; // Delta x
        if ( xmin >= xmax || n < 2 ) {
            reject("invalid parameters passed to seq function");
        }
        mesh = (xmax - xmin)/(n-1);
        return xmin - mesh + cumulative_sum(rep_vector(mesh, n));
    }
    // First cumulant of V_t
    real kappa1(real t, real lambda, real v0, real g) {
        return lambda * v0 / g * (exp(g*t) - 1);
    }
    // derivative of kappa_1
    real dkappa1(real t, real lambda, real v0, real g) {
        return lambda * v0 * exp(g*t);
    }
    // second cumulant of V_t
    real kappa2(real t, real lambda, real v0, real g) {
        return lambda * v0*v0 / (2*g) * (exp(2*g*t) - 1);
    }
    // derivative of kappa_2
    real dkappa2(real t, real lambda, real v0, real g) {
        return lambda * v0*v0 * exp(2*g*t);
    }
    // standard deviation of sd(V_t) = sqrt(kappa_2)
    real sqrtkappa2(real t, real lambda, real v0, real g) {
        return v0 * sqrt(lambda / (2*g) * (exp(2*g*t) - 1));
    }
    // derivative of sd(V_t)
    real dsqrtkappa2(real t, real lambda, real v0, real g) {
        return 0.5 * dkappa2(t, lambda, v0, g) / sqrtkappa2(t, lambda, v0, g);
    }
    /* auxiliary function z for normalizing constant Z
     * Z = phi(z) where z = sqrt(2*lambda/g) + V0/v0 * sqrt(2*g/lambda).
     * Also used for RNG below.
     */
    real fpt_trunc_const(real lambda, real v0, real g, real V0) {
        return sqrt(2*lambda/g) + V0/v0 * sqrt(2*g/lambda);
    }
    /* log-PDF for the FPT distribution of V_t.
     * This allows us to write a sampling statement as 
     *     tau ~ fpt(ell, lambda, v0, g, 0);
     */
    real fpt_lpdf(real t, real ell, real lambda, real v0, real g, real V0) {
        // probability of reaching the threshold
        real logZ = normal_lcdf(fpt_trunc_const(lambda, v0, g, V0) | 0, 1); 
        real kap1 = kappa1(t, lambda, v0, g);
        real sqkap2 = sqrtkappa2(t, lambda, v0, g);
        real dkap1 = dkappa1(t, lambda, v0, g);
        real dsqkap2 = dsqrtkappa2(t, lambda, v0, g);
        real y = (ell - (kap1 + V0*exp(g*t))) / sqkap2;
        real dydt = -(dkap1 + V0*g*exp(g*t)) / sqkap2 - y * dsqkap2 / sqkap2;
        // conditioned on reaching the threshold
        return normal_lpdf(y | 0, 1) + log(-dydt) - logZ; 
    }
    /* For completeness, we also define the (log) survival function
     * (complementary CDF) for the FPT distribution
     */
    real fpt_lccdf(real t, real ell, real lambda, real v0, real g, real V0) {
        real logZ = normal_lcdf(fpt_trunc_const(lambda, v0, g, V0) | 0, 1);
        real kap1 = kappa1(t, lambda, v0, g);
        real sqkap2 = sqrtkappa2(t, lambda, v0, g);
        real y = (ell - (kap1 + V0*exp(g*t))) / sqkap2;
        return log_diff_exp(normal_lcdf(y | 0, 1), log1m_exp(logZ)) - logZ;
    }
    /* And we define the (log) cumulative distribution function
     */
    real fpt_lcdf(real t, real ell, real lambda, real v0, real g, real V0) {
        real logZ = normal_lcdf(fpt_trunc_const(lambda, v0, g, V0) | 0, 1);
        real kap1 = kappa1(t, lambda, v0, g);
        real sqkap2 = sqrtkappa2(t, lambda, v0, g);
        real y = (ell - (kap1 + V0*exp(g*t))) / sqkap2;
        return normal_lccdf(y | 0, 1) - logZ;
    }
    /* Auxiliary random number generator (RNG) that draws from a 
     * truncated normal distribution. Used for the RNG of the
     * FPT below.
     * This uses a simple rejection scheme.
     */
    real truncnorm_rng(real a, real b, real mu, real sigma) {
        // rejection sampler
        int n = 0;
        real x = normal_rng(mu, sigma);
        while ( x < a || x > b ) {
            x = normal_rng(mu, sigma);
            n = n + 1;
            if ( n > 1000 ) {
                reject("exceeded number of samples during rejection sampling");
            }
        }
        return x;
    }
    /* RNG for the first passage time distribution of V_t.
     * Used in the generated quantities block to simulate 
     * rebound times.
     */
    real fpt_rng(real ell, real lambda, real v0, real g, real V0) {
        real t;
        if ( V0 >= ell ) {
            // trivial case...
            t = 0;
        } else {
            real u = fpt_trunc_const(lambda, v0, g, V0);
            real y = truncnorm_rng(negative_infinity(), u, 0, 1);
            real a = y^2 * lambda * v0^2 / (2*g);
            real b = ell + lambda * v0 / g;
            real c = lambda * v0 / g + V0;
            real D = a*(a - c^2 + b^2);
            real x; // x solved from root formula for quadratic eqn
            if ( y > 0 ) {
                x = (b*c + sqrt(D))/(c^2-a);
            } else {
                x = (b*c - sqrt(D))/(c^2-a);
            }
            // V0 < ell ensures that the discriminant is positive
            t = log(x)/g;
        }
        return t;
    }
}

data {
    int<lower=0> NumSubjects;
    int<lower=0> NumTimePts[NumSubjects]; // number of VL observations per subject
    vector[max(NumTimePts)] TimePts[NumSubjects]; // observation times for each subject
    int<lower=0> NumSimTimePts; // used in generated quantities for simulation of VL
    vector[max(NumTimePts)] VirusLoad[NumSubjects]; // VL observations
    int CensorCode[NumSubjects, max(NumTimePts)]; // VL censoring
    real<lower=0> StartART[NumSubjects]; // used as covariate for growth rate and lambda
    real<lower=0> DetectionLimit; // the LoD of the VL, required for computing reboundTime
    // hyper parameters for priors
    real PriorMeanLogK; // hyper-parameter for carying capacity K
    real<lower=0> PriorSdLogK; // hyper-parameter for carying capacity K
    real PriorMeanLogR; // hyper-parameter for growth rate r
    real<lower=0> PriorSdLogR; // hyper-parameter for growth rate r
    real<lower=0> PriorMeanAlphaLogR; // hyper-parameter for effect of ART init on growth rate
    real<lower=0> PriorSdAlphaLogR; // hyper-parameter for effect of ART init on growth rate
    real<lower=0> MaxR; // set an upper bound for the growth rate to avoid numerical issues
    real<lower=0> PriorMeanSigma; // hyper-parameter for VL measurement error
    real<lower=0> PriorSdSigma; // hyper-parameter for VL measurement error
    // priors for the reactivation model
    real PriorMeanLogLambda; // hyper-parameter for recrudescence rate lambda 
    real<lower=0> PriorSdLogLambda;  // hyper-parameter for recrudescence rate lambda
    real PriorMeanLogVZero;  // hyper-parameter for initial VL v0
    real<lower=0> PriorSdLogVZero;  // hyper-parameter for initial VL v0
    real MaxLogVZero; // set an upper bound for v0 to avoid numerical issues
    real PriorMeanAlphaLogLambda; // hyper-parameter for effect of ART init on recrudescence rate
    real<lower=0> PriorSdAlphaLogLambda; // hyper-parameter for effect of ART init on recrudescence rate
    // drug washout delay
    real DrugDelay;
}

/* In the transformed data block, we can pre-compute some 
 * functions of the data that are used in the model
 */
transformed data {
    vector[max(NumTimePts)] LogVirusLoad[NumSubjects];
    real LogDetectionLimit = log(DetectionLimit);
    vector[NumSubjects] StartARTStd; // standardized StartART
    real LogMaxR = log(MaxR);
        
    for ( n in 1:NumSubjects ) {
        LogVirusLoad[n][1:NumTimePts[n]] = log(VirusLoad[n][1:NumTimePts[n]]);
    }
    
    StartARTStd = (to_vector(StartART) - mean(StartART)) / sd(StartART);
}

parameters {
    vector<lower=LogDetectionLimit>[NumSubjects] logK;
    vector<upper=LogMaxR>[NumSubjects] logr;
    real<lower=0> sigma;
    
    real mu_logr;
    real<lower=0> sigma_logr;
    real alpha_logr; // weight of StartART
    
    real mu_logK;
    real<lower=0> sigma_logK;
    
    // reactivation model parameters
    real<upper=min({MaxLogVZero, LogDetectionLimit})> logv0;
    vector[NumSubjects] loglambda;
    real mu_loglambda;
    real<lower=0> sigma_loglambda;
    real alpha_loglambda;
    real<lower=0> fstReactivTime[NumSubjects];
    // rebound time relative to first reactivation time
    real<lower=0> relativeReboundTime[NumSubjects];
}

transformed parameters {
    real<lower=0> reboundTime[NumSubjects]; // the actual rebound time
    vector<upper=0>[NumSubjects] logTheta; // log(ell / K)
    vector<upper=0>[NumSubjects] log1mTheta; // log(1-ell/K)
    vector<lower=0>[NumSubjects] r;
    vector<lower=0>[NumSubjects] lambda;
    real<lower=0> v0;
    
    logTheta = LogDetectionLimit - logK;
    for ( n in 1:NumSubjects ) {
        log1mTheta[n] = log1m_exp(logTheta[n]);
    }
    
    r = exp(alpha_logr * StartARTStd + logr);
    lambda = exp(alpha_loglambda * StartARTStd + loglambda);
    
    v0 = exp(logv0);
    
    for ( n in 1:NumSubjects ) {
        reboundTime[n] = DrugDelay + fstReactivTime[n] + relativeReboundTime[n];
    }
}

/* The model block defines the prior density and the likelihood function.
 * Our observations consist of VL measurements.
 * The FPT distribution enters the model as a prior on the (relative)
 * rebound time.
 */
model {
    for ( n in 1:NumSubjects ) {
        logr[n] ~ normal(mu_logr, sigma_logr) T[,LogMaxR];
        logK[n] ~ normal(mu_logK, sigma_logK) T[LogDetectionLimit,];
    }

    sigma ~ normal(PriorMeanSigma, PriorSdSigma);
    
    mu_logr ~ normal(PriorMeanLogR, PriorSdLogR); 
    sigma_logr ~ normal(0.0, PriorSdLogR);
    alpha_logr ~ normal(PriorMeanAlphaLogR, PriorSdAlphaLogR);
    
    mu_logK ~ normal(PriorMeanLogK, PriorSdLogK);
    sigma_logK ~ normal(0.0, PriorSdLogK);
    
    // priors for the reactivation model
    loglambda ~ normal(mu_loglambda, sigma_loglambda);
    mu_loglambda ~ normal(PriorMeanLogLambda, PriorSdLogLambda);
    sigma_loglambda ~ normal(0.0, PriorSdLogLambda);
    alpha_loglambda ~ normal(PriorMeanAlphaLogLambda, PriorSdAlphaLogLambda);
    
    logv0 ~ normal(PriorMeanLogVZero, PriorSdLogVZero);
       
    // likelihood of the VL observations
    for ( n in 1:NumSubjects ) {
        // auxiliary variables
        vector[NumTimePts[n]] ts = TimePts[n][1:NumTimePts[n]];
        vector[NumTimePts[n]] logVLhat = log_logistic_model(ts, r[n], 
                logK[n], logTheta[n], log1mTheta[n], reboundTime[n]);
                
        for ( i in 1:NumTimePts[n] ) {
            LogVirusLoad[n,i] ~ censored_normal(logVLhat[i], sigma, CensorCode[n,i]);
        }
    }
    
    // the first reactivation time is exponentially distributed
    for ( n in 1:NumSubjects ) {
        fstReactivTime[n] ~ exponential(lambda[n]);
    }
    
    // likelihood of the relative rebound time
    for ( n in 1:NumSubjects ) {
        relativeReboundTime[n] ~ fpt(DetectionLimit, lambda[n], v0, r[n], v0);
    }
}

/* The generated quantities block can be used for prediction and simulation
 * We compute the predicted VL and simulate VL observations in order to
 * make figures of the model fit.
 * We also sample from the FPT distribution, and we compute the likelihood of
 * the VL measurements so that we can compute WAIC
 */
generated quantities {
    vector[NumSimTimePts] logVLhat[NumSubjects];
    vector[NumSimTimePts] logVLsim[NumSubjects];
    real timeToReboundsim[NumSubjects];
    real timeToReboundsimCndlFstReactiv[NumSubjects];
    real extrapolatedFstReactivTime[NumSubjects]; // doesn't have to be positive
    vector[sum(NumTimePts)] loglikes;

    for ( n in 1:NumSubjects ) {
        //real logitTheta = logTheta[n] - log1mTheta[n]; // auxiliary parameter
        vector[NumSimTimePts] simTimePts = seq(0, TimePts[n, NumTimePts[n]], NumSimTimePts);
        // simulate VL data, predict VL curves
        vector[NumSimTimePts] predictions = log_logistic_model(simTimePts, r[n], logK[n], 
                logTheta[n], log1mTheta[n], reboundTime[n]);
        logVLhat[n] = predictions;
        for ( i in 1:NumSimTimePts ) {
            logVLsim[n, i] = predictions[i] + normal_rng(0, sigma);
        }
        // sample from the FPT distribution with initial condition V_0 = v0
        timeToReboundsimCndlFstReactiv[n] = DrugDelay + 
            fstReactivTime[n] + 
            fpt_rng(DetectionLimit, lambda[n], v0, r[n], v0);
        /* sample from the FPT distribution with initial condition V_0 = 0
         * by marginalizing over T_1
         */
        timeToReboundsim[n] = DrugDelay + 
            exponential_rng(lambda[n]) +
            fpt_rng(DetectionLimit, lambda[n], v0, r[n], v0);
    }
    
    // compute extrapolated first reactivation time under deterministic exponential growth
    for ( n in 1:NumSubjects ) {
        extrapolatedFstReactivTime[n] = reboundTime[n] - (LogDetectionLimit - logv0)/r[n];
    }
    
    // record log-likelihoods of observations for WAIC computation
    for ( n in 1:NumSubjects ) {
        // auxiliary variables
        vector[NumTimePts[n]] ts = TimePts[n][1:NumTimePts[n]];
        vector[NumTimePts[n]] xs = log_logistic_model(ts, r[n], 
                logK[n], logTheta[n], log1mTheta[n], reboundTime[n]);
        // store loglikes in the loglikes vector
        for ( i in 1:NumTimePts[n] ) {
            loglikes[sum(NumTimePts[:(n-1)])+i] = 
                    censored_normal_lpdf(LogVirusLoad[n,i] | xs[i], sigma, CensorCode[n,i]);
        }
    }
}
