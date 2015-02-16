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

############# 16 feb 2015

#df = pd.read_csv("DRUG_ANALYSIS_SET.csv", index_col=0)


#Look for drug names:
PLX4720 PLX4720 (rescreen) 1036 1371
BMS-345541 203
AZD6244 1062 1498
AZD8055 1059
AZ6482 AZD6482 (wrong written maybe) 156  1066
5Z-7-OXO (5Z)-7-Oxozeaenol 1242
BX-795 1037

drug_name = [x[0] for x in df[['DRUG_NAME']].values]

len(drug_name)
len(set(drug_name))
# Which is the exact drug name?
for name in drug_name:
    if 'plx4720' in name.lower():
        print name

# Which is the index of each chosen drug?
# df[df['DRUG_NAME'].isin(['drug_with_the_exact_name'])]

# 1. PLX4720 / PLX4720 (rescreen)           1036 1371
# 2. BMS-345541                             203
# 3. AZD6244                                1062 1498
# 4. AZD8055                                1059
# 5. AZD6482                                156  1066
# 6.(5Z)-7-Oxozeaenol                       1242
# 7. BX-795                                 1037
#
# df.loc[[1036,1371,203,1062,1498,1059,156,1066,1242,1037]]

#df2 = pd.read_csv("IC50_v17.csv", index_col=0)
# df2.loc[[1036,1371,203,1062,1498,1059,156,1066,1242,1037]]