import pandas as pd
import numpy as np
import scipy.linalg as la
import seaborn as sbn


B = np.matrix([-10, -9, 5, 10])
Yprev = np.linspace(-10, 10, num=10)
Y=np.matrix(Yprev)

def ols_inverse(B):
    return np.dot(la.pinv(np.dot(B.T, B)), B.T)


def estimate_X(B, Y, noise_level=0):
    #shape = B.shape``
    #Bprime = B + np.random.randn(shape[0], shape[1])*noise_level
    Bprime = B
    #U,D,V = la.svd(Xprime, full_matrices = False)
    #if D.min()<0.001:
    X_hat = np.dot(ols_inverse(Bprime), Y)
    return X_hat


X_hat = estimate_X(B, Y)

print X_hat
