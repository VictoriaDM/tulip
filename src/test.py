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

df[df['Cell.line.name'].isin([
'CCK-81',
'COLO-320-HSR',
'DIFI',
'HCT-116',
'HT-115',
'HT-29',
'SK-CO-1',
'SNU-C2B',
'SNU-C5',
'SW1116',
'SW1463',
'SW620',
'SW837',
'CaR-1'
])]

1240123         CCK-81            CCK-81    1240123  1710STDY5035353
910569    COLO-320-HSR      COLO-320-HSR     910569  1710STDY5035359
1789883           DIFI              DIFI    1789883   CGP_CCL5369320
905936         HCT-116           HCT-116     905936  1710STDY5035336
907289          HT-115            HT-115     907289  1710STDY5035372
905939           HT-29             HT-29     905939  1710STDY5035371
909718         SK-CO-1           SK-CO-1     909718   CGP_CCL5369288
909740         SNU-C2B           SNU-C2B     909740  1710STDY5035360
1674021         SNU-C5            SNU-C5    1674021   CGP_CCL5369323
909746          SW1116            SW1116     909746  1710STDY5035365
909748          SW1463            SW1463     909748  1710STDY5035364
905962           SW620             SW620     905962  1710STDY5035347
909755           SW837             SW837     909755  1710STDY5035348
924108          CaR-1             CaR-1     924108  CGP_CCL5103703

[1240123, 910569, 1789883, 905936, 907289, 905939, 909718, 909740, 1674021, 909746, 909748, 905962, 909755, 924108] #COSMIC.ID for our cell lines of interest.

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
    if 'drug_in_lower_letters' in name.lower():
        print name

# Which is the index of each chosen drug? 10 indexes (we have repeated drugs)
# df[df['DRUG_NAME'].isin(['drug_with_the_exact_name'])]

# 1. PLX4720 / PLX4720 (rescreen)           1036 1371
# 2. BMS-345541                             203
# 3. AZD6244                                1062 1498
# 4. AZD8055                                1059
# 5. AZD6482                                156  1066
# 6.(5Z)-7-Oxozeaenol                       1242
# 7. BX-795                                 1037
#

drugs_id = [1036, 1371, 203, 1062, 1498, 1059, 156, 1066, 1242, 1037]

cells_id = [1240123, 910569, 1789883, 905936, 907289, 905939, 909718, 909740, 1674021, 909746, 909748, 905962, 909755,
            924108]  #COSMIC_ID for our cell lines of interest. 10 cell lines.

cells_id = ['1240123', '910569', '1789883', '905936', '907289', '905939', '909718', '909740', '1674021', '909746',
            '909748', '905962', '909755', '924108']  #COSMIC_ID for our cell lines of interest.


new_df = df[df['COSMIC_ID'].isin(cells_id)]

new_df['DRUG_ID'].isin(drugs_id)
drugs_we_dont_havE = new_df[~new_df['DRUG_ID'].isin(drugs_id)]



df = pd.read_csv("gdsc_drug_sensitivity_fitted_data_w5.csv", index_col=0, low_memory=False)


#Looking for drugs and cells inside IC50_v17.csv:

df = pd.read_csv("IC50_v17.csv", index_col=0)
drugs_id = [1036, 1371, 203, 1062, 1498, 1059, 156, 1066, 1242, 1037]

cells_id = ['1240123', '910569', '1789883', '905936', '907289', '905939', '909718', '909740', '1674021', '909746',
            '909748', '905962', '909755', '924108']  #COSMIC_ID for our cell lines of interest.

cell_filter = df[cells_id] # cell_id=columns.
cell_drug_filter = cell_filter.ix[drugs_id] # drug_id=indexes.

cd.shape

# Number of NaNs: 18
cell_drug_filter.isnull().sum().sum()


# Delete repeated drug names:

# 1. PLX4720 / PLX4720 (rescreen)           1036 1371
# 3. AZD6244                                1062 1498
# 5. AZD6482                                156  1066

drugs_id = [1371, 203, 1498, 1059, 1066, 1242, 1037]

#Drugs with NaNs: 203, 1059, 1066, 1037.
# nothing to do with 1066 and 1037 (as they have already repeated measurements)
# There are no repeated measurements for BMS-345541 and AZD8055.

df = pd.read_csv("IC50_v17.csv", index_col=0)
drugs_id = [1371, 203, 1498, 1059, 1066, 1242, 1037]
cells_id = ['1240123', '910569', '1789883', '905936', '907289', '905939', '909718', '909740', '1674021', '909746',
            '909748', '905962', '909755', '924108']  #COSMIC_ID for our cell lines of interest.
cell_filter = df[cells_id] # cell_id=columns.
cell_drug_filter = cell_filter.ix[drugs_id] # drug_id=indexes.
cell_drug_filter

##             20 February
##