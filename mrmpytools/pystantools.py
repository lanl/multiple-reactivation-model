"""
stan_tools.py
Some functions to augment pystan.
"""

import pystan
import pickle
import os
from hashlib import md5
import numpy as np
import scipy.stats as sts
import warnings
from scipy.optimize import minimize


def cachedStanModel(model_code, model_name=None, recompile=False, work_dir="/tmp/", **kwargs):
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = os.path.join(work_dir, 'cached-model-{}.pkl'.format(code_hash))
    else:
        cache_fn = os.path.join(work_dir, 'cached-{}-{}.pkl'.format(model_name, code_hash))
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        recompile = True
    if recompile:
        if model_name is None:
            sm = pystan.StanModel(model_code=model_code, **kwargs)
        else:
            sm = pystan.StanModel(model_code=model_code, model_name=model_name, **kwargs)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        print("Using cached StanModel")
    return sm


def log_sum_exp(xs):
    u = np.max(xs)
    return u + np.log(np.sum(np.exp(np.array(xs)-u)))


def log_mean_exp(xs):
    """simply log_sum_exp - log(N) where N is the length of the vector xs"""
    return log_sum_exp(xs) - np.log(len(xs))


def calcWAIC(log_like):
    """
    WAIC = \widehat{lpd} - p_{WAIC} with
    \widehat{lpd} = \sum_n \log(1/S \sum_s p(D_n|theta_s)) and
    p_{WAIC} = \sum_n V_s[\log(p(D_n|theta_s))]
    where V is the sample variance (i.e. with a factor 1/(S-1))
    NB: the input is supposed to be of shape (N, S), where
    N is the number of observations, and S is the number
    of posterior samples.
    """
    N, S = len(log_like), len(log_like[0])
    if S < N:
        warnings.warn("The number of samples is smaller than the number of observations. Perhaps the input matrix should have been transposed.")
    ## compute stuff...
    hat_lpd = np.sum([log_mean_exp(ll) for ll in log_like])
    p_waic = np.sum([np.var(ll, ddof=1) for ll in log_like])
    waic = -2 * (hat_lpd - p_waic)
    return waic, hat_lpd, p_waic


def approximate_mode(xs):
    """
    use a gaussian kernel density to find the location of the mode from
    a discrete sample of a distribution
    """
    kde = sts.gaussian_kde(xs)
    res = minimize(lambda x: -kde(x), np.mean(xs), bounds=[(np.min(xs), np.max(xs))])
    return res.x[0]
