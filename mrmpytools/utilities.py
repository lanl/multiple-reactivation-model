"""
utilities.py
Some common utility functions.
"""

def unzip(xs):
    if len(xs) == 0:
        raise Exception("can't unzip an empty list")
    n = len(xs[0])
    return [[x[i] for x in xs] for i in range(n)]

def unique(xs):
    return sorted(list(set(xs)))