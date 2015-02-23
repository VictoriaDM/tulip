"""This script reads a MIDAS from colon data (CORDA)
replace the data with zscores. Extracts only data 
with inhibitors and no stimuli, drops time 0 and 
save the final data into a json file.

It creatse the 14 json files that are provided in ./data
and can be used as input to the regression algos

"""
# needs CORDA or cno XMIDAS function.


from corda import pl2, cellLines

p = pl2.CORDAPreprocessing(path='../data')
p.fold_change_to_pvalue()


for cell in cellLines:
    print(cell)
    # read a MIDAS file for the structure only
    m = XMIDAS("MD-%s_CORDAnormed_1.csv" % cell)

    # replace data with pvalues (actually zscores)
    m.df = p.pvalues[cell].copy()
    indices = m.stimuli.index[m.stimuli.sum(axis=1) == 0]

    # remove experiments keeping only those with no stimuli  but 1 inihibitor
    m.remove_experiments([x for x in m.experiments.index if x not in indices])

    # no need for time 0 (FC=0)
    m.remove_times(0)

    # let us drop the indices (cell/exp/time)
    df = m.df.copy()
    df.reset_index(drop=True, inplace=True)

    # and rename the index basde on the inhibitor
    df.index = ['PI3K', 'mTOR', 'MEKi', 'TBK1', 'IKK', 'PLX', 'TAK']

    # finally, we save the file
    df.to_json('%s_zscores.json' % cell)


