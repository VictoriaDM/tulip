import pandas as pd
import numpy as np
import scipy.linalg as la

# CREATING B and Y.
B = np.matrix([-10, -9, 5, 10])    # B = 4 values.
Y = np.matrix(np.linspace(-10, 10, num=20))  # Y = 20 IC50s.

def ols_inverse(B):
    return np.dot(la.pinv(np.dot(B.T, B)), B.T)

def estimate_X(B, Y, noise_level=0):
    shape = B.shape
    Bprime = B + np.random.randn(shape[0], shape[1])*noise_level
    #U,D,V = la.svd(Xprime, full_matrices = False)
    #if D.min()<0.001:
    X_hat = np.dot(ols_inverse(Bprime), Y)
    return X_hat

X_hat = estimate_X(B, Y)

print X_hat

# ! ONLY CHANGE = B is transposed !
B = B.T
X_hat = X_hat.T

p = X_hat.shape[1]   # len(X.columns)
N = X_hat.shape[0]   # len(df)


# PREDICTING Y using B and X_hat (calculated using ols_inverse and estimate_X functions)
def linear_model(X, beta):
    """Returns Y = beta . X
    :param X: a matrix of shape (N,p)
    :param list beta: list of p-length

    """
    return np.dot(X, beta)

def predict(X_hat, B, noise_level=0):
    this = linear_model(X_hat+np.random.randn(N,p)*noise_level, B)
    return this

yerrors = []
for i in range(0,10000):
    y = predict(X_hat, B, noise_level=0.0)
    yerrors.append(y)
yerrors = pd.DataFrame(yerrors)

# Question: Why I cannot convert yerrors into a dataframe?
# Python output Error: "ValueError: Must pass 2-d input"
