"""
Some notebooks/literature on regression
===========================================

* http://nbviewer.ipython.org/github/temporaer/tutorial_ml_gkbionics/blob/master/3b%20-%20Linear%20regression%202D.ipynb
* http://nbviewer.ipython.org/github/jbpoline/bayfmri/blob/master/notebooks/005-Simple-Linear-Regression.ipynb
"""

import numpy as np


def linear_model(X, beta):
    """


    :param X: a matrix of shape (N,p)
    :param list beta: list of p-length

    """
    return np.dot(X, beta)


def ols(X):
    """A simple implementation of Ordinary Least Square method

    :param X: matrix which solves the OLS regression problem for full-rank X
    """
    return np.dot(la.pinv(np.dot(X.T, X)), X.T)






class Regression(object):
    pass
