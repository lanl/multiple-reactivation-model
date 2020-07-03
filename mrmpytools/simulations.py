"""
simulations.py
Simulate the multiple-reactivation model.
"""

import scipy.stats as sts
import numpy as np
import scipy.optimize
import scipy.special


def sampleVpath(tmax, g, v0, lam0):
    """
    Sample trajectories of the V_t process.
    
    Args:
        tmax: simulate until this time
        g: the exponential growth rate
        v0: the initial viral load after a successful reactivation event
        lam0: the recrudescence rate
        
    Returns:
        Vtots: tuples of jump times and total viral load at these times.
            both with left and right limit of the jump process to 
            facilitate plotting
        Ts: the jump times. Convenient for plotting individual clones.
    """
    tn = 0.0
    Is = []
    Ts = []
    while tn < tmax:
        dt = sts.expon.rvs(scale = 1/lam0)
        if tn + dt < tmax:
            Is.append(dt)
        else:
            Is.append(tmax-tn)
        Ts.append(np.sum(Is))
        tn = Ts[-1]        
    Vss = [[v0 * np.exp(g*(t-T)) for T in Ts if t >= T] for t in Ts]
    lVtots = [(t, np.sum(Vs[:-1])) for t, Vs in zip(Ts, Vss) if len(Vs) > 0]
    rVtots = [(t, np.sum(Vs)) for t, Vs in zip(Ts, Vss) if len(Vs) > 0]
    Vtots = [val for tup in zip(lVtots, rVtots) for val in tup]
    return Vtots, Ts



def sampleFPT(ell, g, v0, lam0, nmax=100):
    """
    Sample an exact rebound time using the process V_t.
    The rebound time is the first passage time of the 
    limit of detection.
    
    Args:
        ell: limit of detection
        g: exponential growth rate
        v0: initial viral load after successful reactivation
        lam0: recrudescence rate
    
    Kwargs:
        nmax: max number of reactivations (for safety)
    """
    t = 0
    V = 0;
    for n in range(nmax):
        ## sample the time of the next event
        t_next = t + sts.expon.rvs(scale=1/lam0)
        V_next_minus = V * np.exp(g*(t_next-t)) ## limit from the left
        V_next = V_next_minus + v0 ## limit from the right
        if V_next >= ell:
            ## find the exact crossing time
            if V_next_minus < ell:
                return t_next
            else:
                return t + np.log(ell/V)/g                            
        V = V_next
        t = t_next
    raise Exception("max iter reached")

    
## simulate the MRM with variation in g

def sampleVGpath(tmax, ggen, v0, lam0, n=1000):
    """
    Simulate the generalized multiple-reactivation model
    with variation in the growth rate (G).
    In this case, the trajectory of V_t between reactivation
    events is not a straight line (on the log scale), 
    because each clone has a different growth rate.
    Therefore, we specify the number of sampling points.
    
    Args:
        tmax: simulate until this time
        ggen: a callable object that returns a (random) growth rate
        v0: the initial viral load after a successful reactivation
        lam0: the recrudescence rate
    Kwargs:
        n: the number of sampling points of the trajectory
    Returns:
        Vtots: tuples of the sampling times and the total viral load
        Ts: the recrudescence times
        gs: the sampled growth rates from ggen
    """
    ts = np.linspace(0, tmax, n)    
    tn = 0.0
    Is = []
    Ts = []
    gs = []
    while tn < tmax:
        dt = sts.expon.rvs(scale = 1/lam0)
        if tn + dt < tmax:
            Is.append(dt)
        else:
            Is.append(tmax-tn)
        Ts.append(np.sum(Is))
        gs.append(ggen())
        tn = Ts[-1]
    Vss = [[v0 * np.exp(g*(t-T)) for g, T in zip(gs, Ts) if t >= T] for t in ts]
    Vtots = [(t, np.sum(Vs)) for t, Vs in zip(ts, Vss) if len(Vs) > 0]
    return Vtots, Ts, gs


def sampleGReboundTime(ell, ggen, v0, lam0):
    """
    Sample first passage time of the limit of detection of the process V_t
    with variable growth rate (G). As the trajectories are not
    exponential between events (as in the case of a fixed g),
    we have to find the time that V_t reaches the LoD numerically.
    We use the root finding algorithm from scipy.optimize
    with the function log_sum_exp(x)-ell.
    
    Args:
        ell: the limit of detection
        ggen: a callable object returning a (random) growth rate
        v0: the initial viral load after a successful reactivation
        lam0: the recrudescence rate
    Returns:
        a randomly sampled reboud time.
    """
    Is = []
    Ts = []
    gs = []
    while True:
        dt = sts.expon.rvs(scale = 1/lam0)
        ## compute the time V reaches ell if no addl event occurs
        if len(Ts) > 0:
            fun = lambda t: scipy.special.logsumexp([g*(t-T) for g, T in zip(gs, Ts)]) - np.log(ell/v0)
            sol = scipy.optimize.root(fun, Ts[-1])
            tau = sol.x[0]
            if tau < Ts[-1] + dt:
                return tau
        ## else...
        Is.append(dt)
        Ts.append(np.sum(Is))
        gs.append(ggen())
        if scipy.special.logsumexp([0] + [g*(Ts[-1]-T) for g, T in zip(gs, Ts)]) > np.log(ell/v0):
            return Ts[-1]    