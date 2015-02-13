import pandas as pd

# cd /Users/Sandman/Desktop/Drugsfeb2015
df = pd.read_csv("IC50_v17.csv")
df2 = pd.read_csv("DRUG_ANALYSIS_SET.csv")
df3 = pd.read_csv("MASTER_LIST.csv")

#df[['Cell.line.name']].transpose().values


## LOOKING CORRECT NAMES FROM THOMAS LIST:

#df = pd.read_csv("MASTER_LIST.csv", index_col=0)
#cell_lines = [x[0] for x in df[['Cell.line.name']].values]

## To check only unique values:
#len(cell_lines)
#len(set(cell_lines))

#for name in cell_lines:
 #   if 'nameofcelllineofinterest' in name.lower():
  #      print name

#CORRECT CELL LINE NAMES:
#CCK-81
#COLO-320-HSR
#DIFI
#HCT-116
#HT-115
#HT-29
#SK-CO-1
#SNU-C2B
#SNU-C5
#SW1116
#SW1463
#SW620
#SW837
