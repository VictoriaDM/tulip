
cellLines = ['CCK81', 'COLO320HSR', 'DIFI', 'HCT116', 'HT115', 'HT29', 'SKCO1', 'SNUC2B', 'SNUC5', 'SW1116', 'SW1463', 'SW620', 'SW837'] #13 cell lines?



data = {cell: pd.read_json(open('%s_zscores.json' % cell,'r').read()) for cell in cellLines}

#cellrange = range(14)
df = pd.DataFrame(columns=cellLines, index=readouts)
inhibitors= ['PI3K', 'mTOR', 'MEKi', 'TBK1', 'IKK', 'PLX', 'TAK']
readouts= ['AKT', 'EGFR', 'ERK', 'GSK3', 'IkBa', 'JNK', 'MEK', 'PI3Kp85', 'RPS6', 'RSKp90', 'SMAD2', 'cJun', 'mTOR', 'p38' ]

for inh in inhibitors:
    for cell in cellLines:
        xval = data[cell]
        df[cell] = xval.ix[inh]
    df.to_json('%s_zscores_inh.json' % inh)


 # CHECK:


new_inh = pd.read_json(open('TAK_zscores_inh.json', 'r').read())
            CCK81  COLO320HSR       DIFI    HCT116     HT115      HT29  \
AKT     -9.186348   -0.138670 -10.340835 -5.730009 -1.984847 -0.369769
EGFR    -0.132675   -0.291535  -3.476061 -1.272318  1.502578  0.273333
ERK     -0.161238   -0.521910  -2.001923 -4.352952  0.477095  0.312730
GSK3    -1.052448   -1.566393  -1.617044 -2.270345 -1.428935 -0.475313

old_cell = pd.read_json(open('DIFI_zscores.json', 'r').read())

AKT       EGFR       ERK      GSK3      IkBa       JNK       MEK  \

TAK  -10.340835  -3.476061 -2.001923 -1.617044 -1.717445 -1.894770  2.191798