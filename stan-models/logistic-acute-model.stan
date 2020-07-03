/* Logistic growth model to estimate 
 * growth rate during acute infection
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
    vector log_logistic_model(vector ts, real r, real logK, 
            real logTheta, real log1mTheta, real tm) {
        vector[num_elements(ts)] xs;
        for ( j in 1:num_elements(ts) ) {
            xs[j] =  logTheta + logK - log_sum_exp(-r*(ts[j]-tm) + log1mTheta, logTheta);
        }
        return xs;
    }
    vector seq(real xmin, real xmax, int n) { // n should be > 2
        real mesh; // Delta x
        if ( xmin >= xmax || n < 2 ) {
            reject("invalid parameters passed to seq function");
        }
        mesh = (xmax - xmin)/(n-1);
        return xmin - mesh + cumulative_sum(rep_vector(mesh, n));
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
    real<lower=0> StartART[NumSubjects];
    real<lower=0> DetectionLimit;
    real PriorMeanLogK;
    real<lower=0> PriorSdLogK;
    real PriorMeanLogR;
    real<lower=0> PriorSdLogR;
    real<lower=0> MaxR;
    real<lower=0> PriorMeanSigma;
    real<lower=0> PriorSdSigma;
    real<lower=0> PriorMeanTau;
    real<lower=0> PriorSdTau;
}

transformed data {
    // Transformed Data
    vector[max(NumTimePts)] LogVirusLoad[NumSubjects];
    real LogDetectionLimit = log(DetectionLimit);
    real LogMaxR = log(MaxR);
        
    for ( n in 1:NumSubjects ) {
        LogVirusLoad[n][1:NumTimePts[n]] = log(VirusLoad[n][1:NumTimePts[n]]);
    }
}

parameters {
    vector<lower=LogDetectionLimit>[NumSubjects] logK;
    vector<upper=LogMaxR>[NumSubjects] logr;
    real<lower=0> sigma;
    
    real mu_logr;
    real<lower=0> sigma_logr;
    
    real mu_logK;
    real<lower=0> sigma_logK;
    
    // the time VL first crosses the LoD
    real<lower=0> tau[NumSubjects];
    real<lower=0> mu_tau;
    real<lower=0> sigma_tau;
}

transformed parameters {
    vector<upper=0>[NumSubjects] logTheta; // log(ell / K)
    vector<upper=0>[NumSubjects] log1mTheta; // log(1-ell/K)
    vector<lower=0>[NumSubjects] r;
    
    logTheta = LogDetectionLimit - logK;
    for ( n in 1:NumSubjects ) {
        log1mTheta[n] = log1m_exp(logTheta[n]);
    }
    
    r = exp(logr);
}

model {
    for ( n in 1:NumSubjects ) {
        logr[n] ~ normal(mu_logr, sigma_logr) T[,LogMaxR];
        logK[n] ~ normal(mu_logK, sigma_logK) T[LogDetectionLimit,];
        tau[n] ~ normal(mu_tau, sigma_tau) T[0,];
    }

    sigma ~ normal(PriorMeanSigma, PriorSdSigma);
    
    mu_logr ~ normal(PriorMeanLogR, PriorSdLogR); 
    sigma_logr ~ normal(0.0, PriorSdLogR);
    
    mu_logK ~ normal(PriorMeanLogK, PriorSdLogK);
    sigma_logK ~ normal(0.0, PriorSdLogK);
    
    mu_tau ~ normal(PriorMeanTau, PriorSdTau);
    sigma_tau ~ normal(0.0, PriorSdTau);
           
    // likelihood of the VL observations
    for ( n in 1:NumSubjects ) {
        // auxiliary variables
        vector[NumTimePts[n]] ts = TimePts[n][1:NumTimePts[n]];
        vector[NumTimePts[n]] logVLhat = log_logistic_model(ts, r[n], 
                logK[n], logTheta[n], log1mTheta[n], tau[n]);
                
        for ( i in 1:NumTimePts[n] ) {
            LogVirusLoad[n,i] ~ censored_normal(logVLhat[i], sigma, CensorCode[n,i]);
        }
    }
}

generated quantities {
    vector[NumSimTimePts] logVLhat[NumSubjects];
    vector[NumSimTimePts] logVLsim[NumSubjects];
    vector[sum(NumTimePts)] loglikes;

    for ( n in 1:NumSubjects ) {
        vector[NumSimTimePts] simTimePts = seq(0, TimePts[n, NumTimePts[n]], NumSimTimePts);
        // simulate VL data, predict VL curves
        vector[NumSimTimePts] predictions = log_logistic_model(simTimePts, r[n], logK[n], 
                logTheta[n], log1mTheta[n], tau[n]);
        logVLhat[n] = predictions;
        for ( i in 1:NumSimTimePts ) {
            logVLsim[n, i] = predictions[i] + normal_rng(0, sigma);
        }
    }
        
    // record log-likelihoods of observations for WAIC computation
    for ( n in 1:NumSubjects ) {
        // auxiliary variables
        vector[NumTimePts[n]] ts = TimePts[n][1:NumTimePts[n]];
        vector[NumTimePts[n]] xs = log_logistic_model(ts, r[n], 
                logK[n], logTheta[n], log1mTheta[n], tau[n]);
        // store loglikes in the loglikes vector
        for ( i in 1:NumTimePts[n] ) {
            loglikes[sum(NumTimePts[:(n-1)])+i] = 
                    censored_normal_lpdf(LogVirusLoad[n,i] | xs[i], sigma, CensorCode[n,i]);
        }
    }
}