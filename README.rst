tulip
#######

todo list
============


* [x] Get data, names of the relevant files
* [x] Use Pandas to read the files and instropect the data
* [x] Plot the IC50 (histogram, 2D plots ...): no need only 1 value used
* [x] identifiy the cell lines 
* [x] identify identifiers of the cell lines
* [x] identifiers of the drug
* [x] identifiy identifers of the targets of interests
* [x] do we have other ICYY ? maybe but private
* [x] write regression algos to identify beta parameters of a linear model
* [x] prototyping with notebook of estimation of beta variation and Y variation in presence of noise on X
* [x] Study of variation of Y and beta as a function of (1) noise level (2) cell lines
* [x] need all JSON of zscores (TC)
* [ ] Look at other regression method
* [ ] Table of Rsquare for the 7 inhibitors for different method e.g. OLS, lasso
* [x] computeing X given Y and beta keeping first column a unit vector
* [ ] New X ?
* [ ] using simulation, compute R2 as functionn of sigma. Where do real data stands ?


* [ ] Cross validation to estimate the mean errors of the test data. Parameters are
    * feature selection: we cannot cover all cases so for N features, we select at most 14 different cases.
      For instance we can select 1 feature out of 14 (14 combi), 2 features (91 combo) and so on but each time, we take only 14 cases.
    * we use 12 cell lines out of 14 for trainign
    * we need to repeat this for the 7 inhibitors
    * we will use 1000 iterations 
    * We will use Linear Regression and Lasso methods
    * This corresponds to 14 times 14 times 7 times 2 methods = 2744 simulations.
    
    We should get a 
    * 2 matrices (for each method) of mean errors of shape 14*14 rows times 7 inhibitors
    * 2 matrices (for each method) of R2 of shape 14*14 rows times 7 inhibitors
      






Notes
======

cell lines
--------------
see src/metadata.py file

Targets/inhibitors
------------------------

Stimuli 
-----------
