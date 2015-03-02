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
            self.IC50 = 1 - pd.read_json(os.sep.join([self.path, os.sep, 'AUC.json']))
        else:
            self.IC50 = pd.read_json(os.sep.join([self.path, os.sep, 'IC50.json']))

        self.inhibitors = self.IC50.index

        self._inhibitor = inhibitor
        self._Y = None
        self._X = None
        self.fillna = 0.84
        self.inhibitor = inhibitor # this set the data and y
        self.method = 'ridge'
        self.l1_ratio = 0.5


    def _set_y(self):
        self._Y =  self.IC50.ix[self.inhibitor].copy()
        self._Y.fillna(self.fillna, inplace=True)
        self._Y = self.Y.values

    def _set_data(self, inhibitor=None):
        if self.inhibitor == 'MEK':
            inh = 'MEKi'
        else:

            inh=self.inhibitor
        filename = os.sep.join([self.path, os.sep,
            'inhibitors/%s_zscores_inh.json' % inh])
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
            progress=True, norm=False, randomise_X=False):
        """For a given inhibitors, and for a given set of features (phosphos)
        return errors for all combo of train/test sets using N=12 cell lines 
        as training and 2 as test.


        """


        from sklearn.linear_model import ElasticNet, Ridge, LinearRegression, Lasso
        if self.method == 'ridge':
            ridge = Ridge(alpha=alpha, normalize=norm)
        elif self.method == 'linear':
            ridge = LinearRegression(normalize=norm)
        elif self.method == 'lasso':
            ridge = Lasso(alpha=alpha, normalize=norm)
        elif self.method == 'elastic':
            ridge = ElasticNet(alpha=alpha, l1_ratio=self.l1_ratio, normalize=norm)
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
        parameters = []
        if progress is True:
            pb = PB(iterations)

        for i in range(0,iterations):
            if random_cells is True:
                Xtrain, Ytrain, Xtest, Ytest, selection = self.train(phosphos, N=12)
            else:
                if randomise_X:
                    X = self.get_randomise_X()
                else:
                    X = self.X.copy()
                X = X[phosphos]
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
            parameters.append(ridge.coef_)
            selections.append([this for this in self.cellLines if this not in
                selection])
            if progress is True:
                pb.animate(i, time.time()-pb.start)

        df = pd.DataFrame({
            'errors':list(flatten(errors)),
            'cellLine': list(flatten(selections))})
        return df, scores, parameters


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


    def get_results(self, alpha=0.01, N=12 , method='ridge', max_combo=14,
            norm=False, randomise_X=False):

        self.method = method
        results = {}
        combo_phosphos = self.get_phospho_combos(M=max_combo)

        combo_phosphos = combo_phosphos[0:]
        pb = PB(len(combo_phosphos)*len(self.inhibitors))

        for j,inhibitor in enumerate(self.inhibitors):
            self.inhibitor = inhibitor

            results[inhibitor] = {'i':[], 'combo':[],'errors':[], 'scores':[], 'params':[]}
            for i, combo in enumerate(combo_phosphos):
                df, scores, params = self.get_errors(list(combo), alpha=alpha, N=N,
                        progress=False, norm=norm, randomise_X=randomise_X)
                results[inhibitor]['errors'].append(df['errors'].mean())
                results[inhibitor]['scores'].append(np.mean(scores))

                pb.animate(j*len(combo_phosphos)+i,time.time() - pb.start)

        df = pd.DataFrame(results)

        scores = pd.DataFrame(dict([(k,df[k]['scores']) for k in df.columns]))
        scores.index = combo_phosphos

        errors = pd.DataFrame(dict([(k,df[k]['errors']) for k in df.columns]))
        errors.index = combo_phosphos

        return scores, errors

    def plot_results(self, scores, errors, save=False, tag=""):
        import pylab
        pylab.figure(1)
        pylab.clf()
        for inh in scores.columns:
            pylab.plot(scores[inh], errors[inh], 'o', 
                    label=inh)
        pylab.legend(loc='best')
        pylab.grid()
        pylab.xlabel("$R^2$ (scores)", fontsize=16)
        pylab.ylabel("Percentage errors", fontsize=16)
        pylab.title("%s regression" % self.method.capitalize())
        if save:
            pylab.savefig("tulip_reg_%s_%s.png" % (self.method, tag))

    def get_randomise_X(self):
        data = self.X.copy()[self.X.columns[1:]]
        for view in np.rollaxis(data.values, 1):
            np.random.shuffle(view)
        return data

    def plot_errors_vs_cv_auc(self, errors, yshift=2, xshift=-1):
        import pylab
        X = self.IC50.std(axis=1)/self.IC50.mean(axis=1)*100
        Y=errors.mean()
        from sklearn.linear_model import LinearRegression
        model = LinearRegression()
        model.fit(X.as_matrix().reshape(7,1), Y.as_matrix())
        print("R2 = ", model.score(X.as_matrix().reshape(7,1), Y))
        pylab.clf()
        pylab.plot(X, Y, 'o')
        pylab.plot(X.values, [x*model.coef_ for x in X.values], 'r')
        pylab.grid()
        pylab.xlabel('CV on AUC for each Drug (across cell lines)', fontsize=16)
        pylab.ylabel('Mean errors (across cell lines)', fontsize=16)
        for i in range(0,7): 
            pylab.text(X[i]+xshift, Y[i]+ yshift, X.index[i])
        pylab.savefig("tulip_errors_function_cv_auc.png")




def get_cv_auc_ic50():
    """AUC is actually 1-AUC"""
    t = Tulip()
    auc = (t.IC50.std()/t.IC50.mean() * 100).mean()
    t = Tulip(use_auc=False)
    ic50 = (t.IC50.std()/t.IC50.mean() * 100).mean()

    return {'cv_ic50': ic50, 'cv_auc':auc}
