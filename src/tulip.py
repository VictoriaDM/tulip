import pandas as pd
import numpy as np
import scipy.linalg as la
#matplotlib.rcParams['figure.figsize'] =(14,8)
from pylab import flatten
import metadata
import linreg
import os
from os.path import join as pj
import itertools
import time

from easydev.progressbar import progress_bar as PB

class Tulip(object):
    def __init__(self, inhibitor='PLX', path='../data/', use_auc=True):
        self.path = path
        if use_auc is True:
            self.IC50 = pd.read_json(os.sep.join([self.path, os.sep, 'AUC.json']))
        else:
            self.IC50 = pd.read_json(os.sep.join([self.path, os.sep, 'IC50.json']))

        self.inhibitors = self.IC50.index

        self._inhibitor = inhibitor
        self._Y = None
        self._X = None
        self.fillna = 0.15
        self.inhibitor = inhibitor # this set the data and y
        self.method = 'ridge'



    def _set_y(self):
        self._Y =  self.IC50.ix[self.inhibitor].copy()
        self._Y.fillna(self.fillna, inplace=True)
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

    def get_errors(self, phosphos, alpha=0.1, N=12, iterations=None,
            progress=True):


        from sklearn.linear_model import Ridge, LinearRegression, Lasso
        if self.method == 'ridge':
            ridge = Ridge(alpha=alpha)
        elif self.method == 'linear':
            ridge = LinearRegression()
        elif self.method == 'lasso':
            ridge = Lasso(alpha=alpha)
        else:
            raise NotImplementedError

        if iterations is None:
            combos = list(itertools.combinations(self.cellLines, N))
            iterations = len(combos)
            random_cells = False
        else:
            random_cells = True


        errors = []
        scores = []
        selections = []
        if progress is True:
            pb = PB(iterations)

        for i in range(0,iterations):
            if random_cells is True:
                Xtrain, Ytrain, Xtest, Ytest, selection = self.train(phosphos, N=12)
            else:
                X = self.X[phosphos]
                selection = list(combos[i])
                indices = [self.cellLines.index(this) for this in combos[i]]

                selection_inv = [this for this in self.cellLines if this not in
                        selection]
                indices_inv = [this for this in range(0,14) if this not in
                        indices]

                Xtrain = X.ix[selection]
                Ytrain = self.Y[indices]
                Xtest = X.ix[selection_inv]
                Ytest = self.Y[indices_inv]


            ridge.fit(Xtrain, Ytrain)
            scores.append(ridge.score(Xtrain, Ytrain))
            percent_error = abs(ridge.predict(Xtest) - Ytest)/abs(Ytest) * 100
            errors.append(percent_error)
            selections.append([this for this in self.cellLines if this not in
                selection])
            if progress is True:
                pb.animate(i, time.time()-pb.start)

        df = pd.DataFrame({
            'errors':list(flatten(errors)),
            'cellLine': list(flatten(selections))})
        return df, scores


    def get_phospho_combos(self, M=14):
        # There are 18383 combos, we select 14 at most for each combo
        # we select i phophos. For each case, we select at most 14 cases
        all_combos = []
        for i in range(1,15):

            combos = list(itertools.combinations(self.phosphos, i))
            if len(combos) > M:
                np.random.shuffle(combos)
                combos = combos[0:M]
            all_combos.extend(combos)
        return all_combos


    def get_results(self, alpha=0.01, N=12 , method='ridge'):

        self.method = method
        results = {}
        combo_phosphos = self.get_phospho_combos()
        pb = PB(len(combo_phosphos)*len(self.inhibitors))
        for j,inhibitor in enumerate(self.inhibitors):
            self.inhbitor = inhibitor
            results[inhibitor] = {'i':[], 'combo':[],'errors':[], 'scores':[]}
            for i, combo in enumerate(combo_phosphos):
                df, scores = self.get_errors(list(combo), alpha=alpha, N=N,
                        progress=False)
                results[inhibitor]['i'].append(i)
                results[inhibitor]['combo'].append(combo)
                results[inhibitor]['errors'].append(df['errors'].mean())
                results[inhibitor]['scores'].append(np.mean(scores))

                pb.animate(j*len(combo_phosphos)+i,time.time() - pb.start)

        df = pd.DataFrame(results)
        return df











