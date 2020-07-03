"""
stats.py
Statistical functions for the multi-reactivation model.
"""

import scipy.stats as sts
import numpy as np
import scipy
from scipy.misc import derivative
import scipy.optimize


def fpt_rng(ell, lam, v0, g, V0):
    """
    random rebound-time generator according to 
    the diffusion approximation of the rebound time
    distribution.
    """
    if V0 >= ell:
        return 0
    ## else...
    u = np.sqrt(2*lam/g) + V0/v0 * np.sqrt(2*g/lam)
    y = sts.truncnorm.rvs(-np.inf, u)
    a = y**2 * lam * v0**2 / (2*g)
    b = ell + lam * v0 / g
    c = lam * v0 / g + V0
    D = a*(a - c**2 + b**2)
    if y > 0:
        x = (b*c + np.sqrt(D)) / (c**2-a)
    else:
        x = (b*c - np.sqrt(D)) / (c**2-a)
    return np.log(x)/g
    
    
def fpt_pdf(t, ell, lam, v0, g, V0):
    """
    Probability density function of the diffusion approximation
    of the rebound time distribution.
    """
    Z = sts.norm.cdf(np.sqrt(2*lam/g) + V0/v0 * np.sqrt(2*g/lam))
    x = np.exp(g*t)
    kap1 = lam*v0/g*(x-1)
    sqkap2 = v0*np.sqrt(lam/(2*g)*(x**2-1))
    dkap1 = lam*v0*x
    dkap2 = lam*(v0*x)**2
    dsqkap2 = 0.5*dkap2/sqkap2
    y = (ell - (kap1 + V0*x)) / sqkap2
    dydt = -(dkap1 + V0*g*x) / sqkap2 - y * dsqkap2 / sqkap2
    return -sts.norm.pdf(y) * dydt / Z; 


## Expectation, standard deviation of the viral load process V_t, 
## and the Pinkevych approxiation

def expectedVt(t, g, v0, lam0):
    """
    Returns the expected viral load E[V_t]
    
    Args:
        t: time post treatment interruption
        g: exponential growth rate
        v0: initial viral load after 
            successful reactivation event
        lam0: the recrudescence rate
    """
    return lam0 * v0 / g * (np.exp(g*t) - 1)



def tildeVt(t, g, v0, lam0, t0):
    """
    Returns the Pinkevych approximation of the viral load V_t
    
    Args:
        t: time post treatment interruption
        g: exponential growth rate
        v0: initial viral load after 
            successful reactivation event
        lam0: the recrudescence rate
        t0: time of first recrudescence event
    """
    return v0*(np.exp(g*(t-t0)) - np.exp(-g/lam0)) / (1-np.exp(-g/lam0))


def sdVt(t, g, v0, lam0):
    """
    Returns the standard deviation of the viral load sqrt(Var[V_t])
    
    Args:
        t: time post treatment interruption
        g: exponential growth rate
        v0: initial viral load after 
            successful reactivation event
        lam0: the recrudescence rate
    """
    return v0 * np.sqrt(lam0 / (2*g) * (np.exp(2*g*t) - 1))



## Fokker-Planck approximation for the FPT distribution


def sigmaV(t, v0, lam, g):
    """
    Standard devation of V_t, same as sdVt
    """
    return v0 * np.sqrt(lam/(2*g) * (np.exp(2*g*t) - 1))

def dsigmaV(t, v0, lam, g):
    """
    Derivative of the standard deviation of V_t
    """
    return 0.5 * v0 * lam * np.exp(2*g*t) / np.sqrt(lam/(2*g) * (np.exp(2*g*t) - 1))

def muV(t, v0, lam, g):
    """
    mean of V_t, same as expectedVt
    """
    return lam * v0 / g * (np.exp(g*t) - 1)

def dmuV(t, v0, lam, g):
    """
    Derivative of the mean of V_t
    """
    return lam * v0 * np.exp(g*t)

def Z(lam, g):
    """
    Normalizing constant Z
    """
    return 1 - sts.norm.cdf(-np.sqrt(2*lam/g))

def FPTpdf(t, ell, v0, lam, g):
    """
    PDF of first passage time distribution of the limit of detection
    """
    sig = sigmaV(t, v0, lam, g)
    mu = muV(t, v0, lam, g)
    dsig = dsigmaV(t, v0, lam, g)
    dmu = dmuV(t, v0, lam, g)
    dydt = -dmu/sig - (ell-mu)/sig * dsig/sig
    return -sts.norm.pdf((ell-mu)/sig) * dydt / Z(lam, g)

def FPTlogpdf(t, ell, v0, lam, g):
    """
    log PDF of first passage time distribution of the limit of detection
    """
    sig = sigmaV(t, v0, lam, g)
    mu = muV(t, v0, lam, g)
    dsig = dsigmaV(t, v0, lam, g)
    dmu = dmuV(t, v0, lam, g)
    dydt = -dmu/sig - (ell-mu)/sig * dsig/sig
    return sts.norm.logpdf((ell-mu)/sig) + np.log(-dydt) - np.log(Z(lam, g))

def FPTcdf(t, ell, v0, lam, g):
    """
    CDF of first passage time distribution of the limit of detection
    """
    sig = sigmaV(t, v0, lam, g)
    mu = muV(t, v0, lam, g)
    z = (ell-mu) / sig
    return (1 - sts.norm.cdf(z)) / Z(lam, g)


def gammaFPTccdf(t, ell, v0, lam, g):
    """
    Complementary CDF (or survival function) for the Gamma-law approximation
    of the first-passage time distribution of the limit of detection.
    """
    a = 2*lam / g *np.tanh(g*t/2)
    x = 2*ell / v0 / (np.exp(g*t)+1)
    return scipy.special.gammainc(a, x)

def gammaFPTpdf(t, ell, v0, lam, g):
    """
    PDF for the Gamma-law approximation of the first-passage time 
    distribution of the limit of detection.
    Computed by numerically differentiating gammaFPTccdf
    """
    return -derivative(gammaFPTccdf, t, args=(ell, v0, lam, g), dx=1e-3)



## functions for the generalized MRM with variable growth rates

def solve_u(sigma, g):
    """
    A convenient distribution for the variable growth rate G
    is given by h -> h/(u * (g-u/2)) on the interval [g-u, g].
    In order to find a width u, given a srandard deviation sigma,
    we have to find a root of an equation.
    
    Args:
        sigma: the target standard deviation of G
        g: the maximum growth rate
    Returns:
        u: the with of the support of G
    """
    fun = lambda u: 12*sigma**2 * (g-u/2)**2 - u**2 * (g**2 - g*u + u**2/6)
    x0 = 2*np.sqrt(3)*sigma
    root = scipy.optimize.root(fun, x0)
    return root.x[0]


def gen_ggen_trunc_propto(sigma, g):
    """
    For simulations, we want to draw random deviates from the distribution
    h -> h / (u * (g-u/2)) on the interval [g-u, g]. We parameterize the
    distribution in terms of the maximum value g, and the standard deviation
    sigma. This function returns a random number generator producing 
    realizations of G.
    
    Args:
        sigma: the standard deviation
        g: the maximum growth rate
    Returns:
        a random number generator for G with maximum g 
        and standard deviation sigma
    """
    u = solve_u(sigma, g)
    def ggen():
        r = sts.uniform.rvs()
        x = np.sqrt(2*u*(g-0.5*u)*r + (g-u)**2)
        return x
    return ggen


def kappaG1(t, g, lam, v0, u):
    """
    First cumulant of the generalized multiple-reactivation model
    """
    if t == 0:
        return 0
    else:
        return v0*lam / (g-u/2) * (np.exp(g*t) * (1-np.exp(-u*t))/(u*t) - 1)
    
def kappaG2(t, g, lam, v0, u):
    """
    Second cumulant of the generalized multiple-reactivation model
    """
    if t == 0:
        return 0
    else:
        return v0**2 * lam / (2*g-u) * (np.exp(2*g*t) * (1-np.exp(-2*u*t))/(2*u*t) - 1)

def gammaFPTccdfG(t, y, v0, lam, g, u):
    """
    Complementary CDF (or survival function) of the rebound time distribution
    for the generalized multiple-reactivation model.
    """
    if t <= 0:
        return 1.0;
    ## else...
    k1 = kappaG1(t, g, lam, v0, u)
    k2 = kappaG2(t, g, lam, v0, u)
    a = k1**2 / k2
    x = y * k1 / k2
    return scipy.special.gammainc(a, x)

def gammaFPTpdfG(t, y, v0, lam, g, u):
    """
    PDF of the rebound time distribution for the generalized 
    multiple-reactivation model. Calculated using numerical
    differentiation of the complementary CDF
    """
    return -derivative(gammaFPTccdfG, t, args=(y, v0, lam, g, u), dx=1e-3)