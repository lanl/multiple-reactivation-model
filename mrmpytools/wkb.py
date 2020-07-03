"""
wkb.py
Functions for computing the WKB approximation of the 
rebound-time distribution.
"""

import numpy as np
import scipy.stats as sts
import scipy.special
import scipy.integrate ## solve_ivp
import scipy.interpolate


def Ein(z):
    """
    Returns Ein(z) = int_0^z (1-exp(-s))/s ds
    
    implemented using other exponential integrals 
    (exp1 or expi, dependent on the sign of z).
    """
    if z == 0:
        return 0.0
    elif z > 0:
        return scipy.special.exp1(z) + np.log(z) + np.euler_gamma
    else:
        return -scipy.special.expi(-z) + np.log(-z) + np.euler_gamma
    
    
def inv_phi_x_lim(t, ell, par):
    """
    Computes p_0 required to arrive at ell at time t.
    Uses Lambert W function.
    """
    g, lam, v0 = par
    tstar = np.log(ell*g / (lam*v0) + 1)/g
    u = lam * v0 / (ell*g) * np.exp(g*t)
    a = -u*np.exp(-u)
    ## k determines the branch of W
    k = 0 if t > tstar else -1
    W = np.real(scipy.special.lambertw(a, k=k, tol=1e-10))
    p0 = -(u + W)
    return p0
    
    
def action_lim(p0, par):
    g, lam, v0 = par
    S = lam * v0 / g * (np.exp(p0) - 1 + Ein(-p0))
    return S
    
    
def normalize(ts, us):
    dts = ts[1:] - ts[:-1]
    ws = 0.5 * (us[1:] + us[:-1])
    return np.dot(dts, ws)


def FPTpdfWKBlim(ts, ell, par):
    g, lam, v0 = par
    p0s = np.array([inv_phi_x_lim(t, ell, par) for t in ts])
    Ss = np.array([action_lim(p0, par) for p0 in p0s])
    Prs = np.array([np.exp(-S/v0) for S in Ss])
    Z = normalize(ts, Prs)
    ZPrs = Prs / Z
    return ZPrs


def FPTpdfWKBinter(ts, fs):
    ## uses the result of FPTpdfWKB
    large_t = ts[-1]*10
    ts_extended = np.pad(ts, (1,1), 'constant', constant_values=(0, large_t))
    fs_extended = np.pad(fs, (1,1), 'constant', constant_values=(0,0))
    f = scipy.interpolate.interp1d(ts_extended, fs_extended, kind='cubic')
    return f

def FPTccdfWKB(f, t):
    ## uses the result of FPTpdfWKB
    res = scipy.integrate.quadrature(f, 0, t, maxiter=1000, tol=1e-10, rtol=1e-10)
    F = 1 - res[0]
    return F
