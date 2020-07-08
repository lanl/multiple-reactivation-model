/* This simplified model is used for sensitivity analyses 
 * with simulated data. Estimate only a single parameter for 
 * the growth rate g (Gsp) i.e. no random effects. idem for lambda
 *
 * Abbreviations:
 * Cndl: conditional
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
    vector log_expon_model(vector ts, real r, real tm, real log_vm) {
        return r*(ts-tm) + log_vm;
    }
    vector seq(real xmin, real xmax, int n) { // n should be > 2
        real mesh; // Delta x
        if ( xmin >= xmax || n < 2 ) {
            reject("invalid parameters passed to seq function");
        }
        mesh = (xmax - xmin)/(n-1);
        return xmin - mesh + cumulative_sum(rep_vector(mesh, n));
    }
    // functions to compute the approximate first-passage time distribution
    real kappa1(real t, real lambda, real v0, real g) {
        return lambda * v0 / g * (exp(g*t) - 1);
    }
    real dkappa1(real t, real lambda, real v0, real g) {
        return lambda * v0 * exp(g*t);
    }
    real kappa2(real t, real lambda, real v0, real g) {
        return lambda * v0*v0 / (2*g) * (exp(2*g*t) - 1);
    }
    real dkappa2(real t, real lambda, real v0, real g) {
        return lambda * v0*v0 * exp(2*g*t);
    }
    real sqrtkappa2(real t, real lambda, real v0, real g) {
        return v0 * sqrt(lambda / (2*g) * (exp(2*g*t) - 1));
    }
    real dsqrtkappa2(real t, real lambda, real v0, real g) {
        return 0.5 * dkappa2(t, lambda, v0, g) / sqrtkappa2(t, lambda, v0, g);
    }
    real fpt_trunc_const(real lambda, real v0, real g, real V0) {
        return sqrt(2*lambda/g) + V0/v0 * sqrt(2*g/lambda);
    }
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
    real fpt_lccdf(real t, real ell, real lambda, real v0, real g, real V0) {
        real logZ = normal_lcdf(sqrt(2*lambda/g) + V0/v0 * sqrt(2*g/lambda) | 0, 1);
        real kap1 = kappa1(t, lambda, v0, g);
        real sqkap2 = sqrtkappa2(t, lambda, v0, g);
        real y = (ell - (kap1 + V0*exp(g*t))) / sqkap2;
        return log_diff_exp(normal_lcdf(y | 0, 1), log1m_exp(logZ)) - logZ;
    }
    real fpt_lcdf(real t, real ell, real lambda, real v0, real g, real V0) {
        real logZ = normal_lcdf(sqrt(2*lambda/g) + V0/v0 * sqrt(2*g/lambda) | 0, 1);
        real kap1 = kappa1(t, lambda, v0, g);
        real sqkap2 = sqrtkappa2(t, lambda, v0, g);
        real y = (ell - (kap1 + V0*exp(g*t))) / sqkap2;
        return normal_lccdf(y | 0, 1) - logZ;
    }
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
    // Data
    int<lower=0> NumSubjects;
    int<lower=0> NumTimePts[NumSubjects];
    vector[max(NumTimePts)] TimePts[NumSubjects];
    int<lower=0> NumSimTimePts;
    vector[max(NumTimePts)] VirusLoad[NumSubjects];
    int CensorCode[NumSubjects, max(NumTimePts)];
    real<lower=0> DetectionLimit; // required for computing reboundTime
    real PriorMeanLogR;
    real<lower=0> PriorSdLogR;
    real<lower=0> MaxR;
    real<lower=0> PriorMeanSigma;
    real<lower=0> PriorSdSigma;
    // priors for the reactivation model
    real PriorMeanLogLambda;
    real<lower=0> PriorSdLogLambda;
    real LogVZero;
}

transformed data {
    // Transformed Data
    vector[max(NumTimePts)] LogVirusLoad[NumSubjects];
    real LogDetectionLimit = log(DetectionLimit);
    real LogMaxR = log(MaxR);
    real VZero = exp(LogVZero);
        
    for ( n in 1:NumSubjects ) {
        LogVirusLoad[n][1:NumTimePts[n]] = log(VirusLoad[n][1:NumTimePts[n]]);
    }
}

parameters {
    real<upper=LogMaxR> logr;
    real<lower=0> sigma;        
    // reactivation model parameters
    real loglambda;
    // rebound time relative to drug washout delay
    real<lower=0> fstReactivTime[NumSubjects];
    // rebound time relative to first reactivation time
    real<lower=0> relativeReboundTime[NumSubjects];
}

transformed parameters {
    real<lower=0> reboundTime[NumSubjects]; // actual rebound time

    real<lower=0> r = exp(logr);
    real<lower=0> lambda = exp(loglambda);
    
    for ( n in 1:NumSubjects ) {
        reboundTime[n] = fstReactivTime[n] + relativeReboundTime[n];
    }
}

model {   
    sigma ~ normal(PriorMeanSigma, PriorSdSigma);
    logr ~ normal(PriorMeanLogR, PriorSdLogR) T[,LogMaxR]; 
    // priors for the reactivation model
    loglambda ~ normal(PriorMeanLogLambda, PriorSdLogLambda);
    
    // likelihood of the VL observations
    for ( n in 1:NumSubjects ) {
        // auxiliary variables
        vector[NumTimePts[n]] ts = TimePts[n][1:NumTimePts[n]];
        vector[NumTimePts[n]] logVLhat = log_expon_model(ts, r, reboundTime[n], LogDetectionLimit);
                
        for ( i in 1:NumTimePts[n] ) {
            // the sampling statement for the VLs
            LogVirusLoad[n,i] ~ censored_normal(logVLhat[i], sigma, CensorCode[n,i]);
        }
    }
        
    // the first reactivation time is exponentially distributed
    for ( n in 1:NumSubjects ) {
        fstReactivTime[n] ~ exponential(lambda);
    }
    
    // likelihood of the relative rebound time
    for ( n in 1:NumSubjects ) {
        relativeReboundTime[n] ~ fpt(DetectionLimit, lambda, VZero, r, VZero);
    }
}

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
        vector[NumSimTimePts] predictions = log_expon_model(simTimePts, r, reboundTime[n], LogDetectionLimit);
        logVLhat[n] = predictions;
        for ( i in 1:NumSimTimePts ) {
            logVLsim[n, i] = predictions[i] + normal_rng(0, sigma);
        }
        // sample from the FPT distribution with initial condition V_0 = v0
        timeToReboundsimCndlFstReactiv[n] = fstReactivTime[n] + 
            fpt_rng(DetectionLimit, lambda, VZero, r, VZero);
        /* sample from the FPT distribution with initial condition V_0 = 0
         * by marginalizing over T_1
         */
        timeToReboundsim[n] = exponential_rng(lambda) +
            fpt_rng(DetectionLimit, lambda, VZero, r, VZero);
    }
    
    // compute extrapolated first reactivation time under deterministic exponential growth
    for ( n in 1:NumSubjects ) {
        extrapolatedFstReactivTime[n] = reboundTime[n] - (LogDetectionLimit - LogVZero)/r;
    }
    
    // record log-likelihoods of observations for WAIC computation
    for ( n in 1:NumSubjects ) {
        // auxiliary variables
        vector[NumTimePts[n]] ts = TimePts[n][1:NumTimePts[n]];
        vector[NumTimePts[n]] xs = log_expon_model(ts, r, reboundTime[n], LogDetectionLimit);
        // store loglikes in the loglikes vector
        for ( i in 1:NumTimePts[n] ) {
            loglikes[sum(NumTimePts[:(n-1)])+i] = 
                    censored_normal_lpdf(LogVirusLoad[n, i] | xs[i], sigma, CensorCode[n, i]);
        }
    }
}
