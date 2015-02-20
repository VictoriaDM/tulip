"""
Some notebooks/literature on regression
===========================================

* http://nbviewer.ipython.org/github/temporaer/tutorial_ml_gkbionics/blob/master/3b%20-%20Linear%20regression%202D.ipynb
* http://nbviewer.ipython.org/github/jbpoline/bayfmri/blob/master/notebooks/005-Simple-Linear-Regression.ipynb
"""

import numpy as np
import scipy.linalg as la


def linear_model(X, beta):
    """Returns Y = beta . X


    :param X: a matrix of shape (N,p)
    :param list beta: list of p-length

    """
    return np.dot(X, beta)


def ols(X):
    """A simple implementation of Ordinary Least Square method

    :param X: matrix which solves the OLS regression problem for full-rank X
    """
    return np.dot(la.pinv(np.dot(X.T, X)), X.T)



def estimate_beta(X, Y, noise_level=0):
    shape = X.shape
    Xprime = X + np.random.randn(shape[0], shape[1])*noise_level
    #U,D,V = la.svd(Xprime, full_matrices = False) 
    #if D.min()<0.001:
    #    print D
    beta_hat = np.dot(ols(Xprime), Y)
    return beta_hat




def ridge_estimate(X, Y, L=10, noise_level=0.1):
    """" 
    Calculate the matrix for ridge regression
    
    Parameters 
    ----------
    X : 2d array
        The design matrix
    
    L : float
        ridge parameter for regularization    
    """
    shape = X.shape
    Xprime = X + np.random.randn(shape[0], shape[1])*noise_level
    Xprime = np.dot(la.pinv(np.dot(Xprime.T,X) + L * np.eye(Xprime.shape[-1])), X.T)
    return np.dot(Xprime, Y)





class Regression(object):
    pass
