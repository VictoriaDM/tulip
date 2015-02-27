import pandas as pd
import numpy as np
import scipy.linalg as la
#matplotlib.rcParams['figure.figsize'] =(14,8)
from pylab import flatten
import metadata
import linreg
import os
from os.path import join as pj


from easydev.progressbar import progress_bar as PB

class Tulip(object):
    def __init__(self, inhibitor='PLX', path='../data/'):
        self.path = path
        self.IC50 = pd.read_json(os.sep.join([self.path, os.sep, 'IC50.json']))
        self.inhibitors = self.IC50.index


        self._inhibitor = inhibitor
        self._Y = None
        self._X = None
        self.inhibitor = inhibitor # this set the data and y
        self.method = 'ridge'



    def _set_y(self,  fillna=1):
        self._Y =  self.IC50.ix[self.inhibitor].copy()
        self._Y.fillna(fillna, inplace=True)
        self._Y = self.Y.values

    def _set_data(self, inhibitor=None):
        filename = os.sep.join([self.path, os.sep,
            'inhibitors/%s_zscores_inh.json' % self.inhibitor])
        df = pd.read_json(open(filename, 'r').read())
        #inhibitor = 'MEK'  # note that for MEK, there is an i in the filename but no in
        # the header...
        df['one'] = [1] * 14
        # put 'one' as first columns
        df = df[['one'] + list(df.columns[0:-1])]
        self.df = df.copy()

    def _get_phosphos(self):
        return [x for x in self.df.columns if x != 'one']
    phosphos = property(_get_phosphos)

    def _get_cellLines(self):
        return [x for x in self.df.index]
    cellLines = property(_get_cellLines)

    def _get_inhibitor(self):
        return self._inhibitor
    def _set_inhibitor(self, inhibitor):
        assert inhibitor in self.inhibitors
        self._inhibitor = inhibitor
        self._set_y()
        self._set_data()
    inhibitor = property(_get_inhibitor, _set_inhibitor)

    def _get_X(self):
        return self.df.copy()
    X = property(_get_X)

    def _get_Y(self):
        return self._Y.copy()
    Y = property(_get_Y)

    def train(self, phosphos, N=10):
        # select rows
        X = self.X[phosphos]

        # select N random cell lines
        selection_integer = range(0, 14)
        np.random.shuffle(selection_integer)
        selection = [self.cellLines[i] for i in selection_integer[0:N]]
        inv = [self.cellLines[i] for i in selection_integer[N:]]

        Xtrain = X.ix[selection]
        Ytrain = self.Y[selection_integer[0:N]]
        Xtest = X.ix[inv]
        Ytest = self.Y[selection_integer[N:]]
        return Xtrain, Ytrain, Xtest, Ytest, selection

    def get_errors(self, phosphos, alpha=0.1, N=12, iterations=100):

        from sklearn.linear_model import Ridge, LinearRegression, Lasso
        if self.method == 'ridge':
            ridge = Ridge(alpha=alpha)
        elif self.method == 'linear':
            ridge = LinearRegression()
        elif self.method == 'lasso':
            ridge = Lasso(alpha=alpha)
        else:
            raise NotImplementedError

        errors = []
        scores = []
        selections = []
        pb = PB(iterations)
        import time
        for i in range(0,iterations):
            Xtrain, Ytrain, Xtest, Ytest, selection = self.train(phosphos, N=12)
            ridge.fit(Xtrain, Ytrain)
            scores.append(ridge.score(Xtrain, Ytrain))
            percent_error = abs(ridge.predict(Xtest) - Ytest)/abs(Ytest) * 100
            errors.append(percent_error)
            selections.append([this for this in self.cellLines if this not in
                selection])
            pb.animate(i, time.time()-pb.start)

        df = pd.DataFrame({
            'errors':list(flatten(errors)),
            'cellLine': list(flatten(selections))})
        return df, scores










