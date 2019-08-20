import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import statistics
from scipy import integrate

files = glob('Y:/experiments/Experiments_004600/004681/After registration/mitotic/local bulk DNA seperation/peakshift_190715_115925_Airyscan_Processed.G2_4681_Mitotic_Hemi_2019_07_15__11_22_34.czi - Image 4 #1.tif_registered.tif/*.csv', recursive=True)
# files = glob('/groups/gerlich/experiments/Experiments_004600/004681/After registration/mitotic/local bulk DNA seperation/peakshift_190715_115925_Airyscan_Processed.G2_4681_Mitotic_Hemi_2019_07_15__11_22_34.czi - Image 4 #1.tif_registered.tif/*.csv', recursive=True)

files


def seperation_bulkDNA_local(dataframe,baseline):
    df = dataframe

    i=0
    if 'Scc1' in df.columns:
        standart = pd.DataFrame()
        standart['Distance'] = df['Distance']
        standart['Scc1'] = (df['Scc1']-df['Scc1'].mean())/df['Scc1'].std()
        standart['f-ara-EdU'] = (df['f-ara-EdU']-df['f-ara-EdU'].mean())/df['f-ara-EdU'].std()
        standart['Hoechst'] = (df['Hoechst']-df['Hoechst'].mean())/df['Hoechst'].std()

        # defiing chromosomal region according to hoechst threshhold
        standart['Chromosome'] = "no"
        # standart['Chromosome_2'] = standart.loc[standart['Hoechst'] > -1] = 'yes'
        standart.loc[standart['Hoechst'] > baseline, ['Chromosome']] = 'yes'

        # splotting values over threshhold into condecutive dictionaries
        s = (standart['Chromosome'] == 'yes')
        s = (s.gt(s.shift(fill_value=False)) + 0).cumsum() * s
        grp = {}
        for i in np.unique(s)[1:]:
            grp[i] = standart.loc[s == i, ['Distance', 'Scc1', 'f-ara-EdU', 'Hoechst']]
    else:
        standart = pd.DataFrame()
        standart['Distance'] = df['Distance']
        # standart['Scc1'] = (df['Scc1']-df['Scc1'].mean())/df['Scc1'].std()
        standart['f-ara-EdU'] = (df['f-ara-EdU']-df['f-ara-EdU'].mean())/df['f-ara-EdU'].std()
        standart['Hoechst'] = (df['Hoechst']-df['Hoechst'].mean())/df['Hoechst'].std()

        # defiing chromosomal region according to hoechst threshhold
        standart['Chromosome'] = "no"
        # standart['Chromosome_2'] = standart.loc[standart['Hoechst'] > -1] = 'yes'
        standart.loc[standart['Hoechst'] > baseline, ['Chromosome']] = 'yes'

        # splotting values over threshhold into condecutive dictionaries
        s = (standart['Chromosome'] == 'yes')
        s = (s.gt(s.shift(fill_value=False)) + 0).cumsum() * s
        grp = {}
        for i in np.unique(s)[1:]:
            grp[i] = standart.loc[s == i, ['Distance', 'f-ara-EdU', 'Hoechst']]

    # extracting the biggest block over threshhold as Chromosomal mass
    dicts = list(range(1, i+1))
    dicts
    names = {}
    for block in dicts:
        names["Block{0}".format(block)] = [block, len(grp[block])]
    currentDF = pd.DataFrame.from_dict(names, orient='index', columns=['grp_number', 'Length'])
    currentDF.head()
    Chromosome_block = currentDF.loc[currentDF['Length'] == currentDF['Length'].max(), 'grp_number']
    # print(Chromosome_block)
    x = int(Chromosome_block)
    chromosome = grp[x]

    # shifting curves above 0
    if (chromosome['f-ara-EdU'].min() < 0):
        chromosome["f-ara-EdU"] = chromosome["f-ara-EdU"] - chromosome['f-ara-EdU'].min()
    if (chromosome['Hoechst'].min() < 0):
        chromosome["Hoechst"] = chromosome["Hoechst"] - chromosome['Hoechst'].min()

    # seperating in left and right sister
    Hoechst_median = statistics.median(chromosome['Distance'])
    left = chromosome.loc[chromosome['Distance'] <= Hoechst_median]
    right = chromosome.loc[chromosome['Distance'] >= Hoechst_median]

    # calculating integrals
    left_Hoechst = integrate.trapz(left['Hoechst'], left['Distance'])
    left_EdU = integrate.trapz(left['f-ara-EdU'], left['Distance'])

    right_Hoechst = integrate.trapz(right['Hoechst'], right['Distance'])
    right_EdU = integrate.trapz(right['f-ara-EdU'], right['Distance'])

    # determining EdU heavy side
    if right_EdU >= left_EdU:
        ratio_EdU = right_EdU/left_EdU
        ratio_Hoechst = right_Hoechst/left_Hoechst
        percentage_EdU = right_EdU/(left_EdU+right_EdU)
        percentage_Hoechst = right_Hoechst/(left_Hoechst+right_Hoechst)
    else:
        ratio_EdU = left_EdU/right_EdU
        ratio_Hoechst = left_Hoechst/right_Hoechst
        percentage_EdU = left_EdU/(left_EdU+right_EdU)
        percentage_Hoechst = left_Hoechst/(left_Hoechst+right_Hoechst)

    results = [left_EdU, right_EdU, ratio_EdU, percentage_EdU, left_Hoechst, right_Hoechst, ratio_Hoechst, percentage_Hoechst]
    return results


def meanpercentage(files,baseline):
    ratios_dict = {}
    index = 0
    for file in files:
        try:
            ratios_dict[f'{index}'] = [file] + seperation_bulkDNA_local(pd.read_csv(file),baseline)
            index += 1
        except KeyError:
            print(f'{file}was skipped for {baseline}')
    ratios_df = pd.DataFrame.from_dict(ratios_dict, orient='index', columns=['File', 'left_EdU', 'right_EdU', 'ratio_EdU', 'percentage_EdU', 'left_Hoechst', 'right_Hoechst', 'ratio_Hoechst', 'percentage_Hoechst'])
    percentage_means = [ratios_df['percentage_EdU'].mean(), ratios_df['percentage_Hoechst'].mean()]
    return percentage_means


baselines = list(np.arange(-1.5, -0.45, 0.05))
print(baselines)
percentage_means_dict = {}

index = 0
for baseline in baselines:
    try:
        percentage_means_dict[index] = [baseline] + meanpercentage(files, baseline)
        index += 1
    except TypeError:
        print(f"{baseline} gave a problem")
percentage_means_df = pd.DataFrame.from_dict(percentage_means_dict, orient='index', columns=['Baseline', 'percentage_EdU', 'percentage_Hoechst'])
percentage_means_df.head


plt.plot(percentage_means_df['Baseline'], percentage_means_df['percentage_EdU'], color='r', label='EdU on one sister')
plt.plot(percentage_means_df['Baseline'], percentage_means_df['percentage_Hoechst'], color='b', label='Hoechst on one sister')
plt.legend()