# import python modules
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import skbio
import session_info

from skbio import DistanceMatrix
from skbio.diversity import beta_diversity, alpha_diversity
from skbio.stats import ordination, distance
from skbio.stats.composition import ancom, clr
from skbio.stats.distance import mantel

from plotnine import ggplot, aes, geom_line, geom_boxplot, geom_jitter, geom_violin
from scipy.stats import mannwhitneyu
from dateutil import rrule

# define functions
# calculate age in weeks
def weeks_between(start_date, end_date):
    weeks = rrule.rrule(rrule.WEEKLY, dtstart=start_date, until=end_date)
    return weeks.count()

rename_ids={'Blank111019': 'NegA', 'Blank15119': 'NegB', 'Blank181019': 'NegC', 'Blank191019': 'NegD', 'Blank201019': 'NegE', 'Blank251019': 'NegF', 'Blank261019': 'NegG', 'Blank271019': 'NegH', 'Blank280919': 'NegI', 'Blank15119': 'NegJ', 'KGCF01': 'X01', 'KGCF02': 'X02', 'KGCF03': 'X03', 'KGCF04': 'X04', 'KGCF05': 'X05', 'KGCF06': 'X06', 'KGCF07': 'X07', 'KGCF10': 'X10', 'KGCF11': 'X11', 'KGCF12': 'X12', 'KGCF13': 'X13', 'KGCF14': 'X14', 'KGCF15': 'X15', 'KGCF16': 'X16', 'KGCF17': 'X17', 'KGCF18': 'X18', 'KGCF19': 'X19', 'KGCF21': 'X21', 'KGCF22': 'X22', 'KGCF24': 'X24', 'KGCF25': 'X25', 'KGCF26': 'X26', 'KGCF27': 'X27', 'KGCF28': 'X28', 'KGCF29': 'X29', 'KGCF30': 'X30', 'KGCF32': 'X32', 'KGCF33': 'X33', 'KGCF34': 'X34', 'KGCF35': 'X35', 'KGCF36': 'X36', 'KGCF37': 'X37', 'KGCF38': 'X38', 'KGCF39': 'X39', 'KGCF40': 'X40', 'KGCF41': 'X41', 'KGCF42': 'X42', 'KGCF43': 'X43', 'KGCF44': 'X44', 'KGCF45': 'X45', 'KGCF46': 'X46', 'KGCF47': 'X47', 'KGCF48': 'X48', 'KGCF49': 'X49', 'KGCF50': 'X50', 'KGCF51': 'X51', 'KGCF52': 'X52', 'KGCF53': 'X53', 'KGCF55': 'X55', 'KGCF56': 'X56', 'KGCF57': 'X57', 'KGCF58': 'X58', 'MCF01s1': 'CF-A-01', 'MCF01s2': 'CF-A-02', 'MCF01s3': 'CF-A-03', 'MCF01s4': 'CF-A-04', 'MCF01s5': 'CF-A-05', 'MCF01s6': 'CF-A-06', 'MCF01s7': 'CF-A-07', 'MCF01s8': 'CF-A-08', 'MCF01s9': 'CF-A-09', 'MCF02s1': 'CF-B-01', 'MCF02s2': 'CF-B-02', 'MCF02s3': 'CF-B-03', 'MCF02s4': 'CF-B-04', 'MCF03s1': 'ZZZ01', 'MCF03s2': 'ZZZ02', 'MCF04s1': 'CF-C-01', 'MCF04s2': 'CF-C-02', 'MCF04s3': 'CF-C-03', 'MCF05s1': 'CF-D-01', 'MCF05s2': 'CF-D-02', 'MCF06s1': 'CF-E-01', 'MCF06s2': 'CF-E-02', 'MCF06s3': 'CF-E-03', 'MCF06s4': 'CF-E-04', 'MCF06s5': 'CF-E-05', 'MCF06s6': 'CF-E-06', 'MCF07s1': 'CF-F-01', 'MCF07s2': 'CF-F-02', 'MCF07s3': 'CF-F-03', 'MCF07s4': 'CF-F-04', 'MCF08s1': 'CF-G-01', 'MCF08s2': 'CF-G-02', 'MCF08s3': 'CF-G-03', 'MCF08s4': 'CF-G-04', 'MCF08s5': 'CF-G-05', 'MCF08s6': 'CF-G-06', 'MCF09s1': 'CF-H-01','MCF09s2': 'CF-H-02', 'MCF09s3': 'CF-H-03', 'MCF09s4': 'CF-H-04', 'MCF09s5': 'CF-H-05', 'MCF09s6': 'CF-H-06', 'MCF10s1': 'CF-I-01', 'MCF10s2': 'CF-I-02', 'MCF10s3': 'CF-I-03', 'MCF10s4': 'CF-I-04', 'MCF11s1': 'CF-J-01', 'MCF11s2': 'CF-J-02', 'MCF11s3': 'CF-J-03', 'MCF11s4': 'CF-J-04', 'MCF12s1': 'CF-K-01', 'MCF12s2': 'CF-K-02', 'MCF12s3': 'CF-K-03', 'MCF12s4': 'CF-K-04', 'MCF12s5': 'CF-K-05', 'MCF13s1': 'CF-L-01', 'MCF13s2': 'CF-L-02', 'MCF13s3': 'CF-L-03', 'MCF13s4': 'CF-L-04', 'MCF14s1': 'CF-M-01', 'MCF14s2': 'CF-M-02', 'WA111019': 'NegK', 'WA121019': 'NegL', 'WA131019': 'NegM', 'WA15119': 'NegN', 'WA161119': 'NegO', 'WA171119': 'NegP', 'WA181019': 'NegQ', 'WA191019': 'NegR', 'WA201019': 'NegS', 'WA271019': 'NegT', 'WA271119': 'NegU', 'WA280919': 'NegV', 'WA281119': 'NegW'}
rename_patient={'MCF01': 'CF-A', 'MCF02': 'CF-B', 'MCF04': 'CF-C', 'MCF05': 'CF-D','MCF06': 'CF-E', 'MCF07': 'CF-F', 'MCF08': 'CF-G', 'MCF09': 'CF-H', 'MCF10': 'CF-I', 'MCF11': 'CF-J', 'MCF12': 'CF-K', 'MCF13': 'CF-L', 'MCF14': 'CF-M', 'KGCF01': 'X', 'KGCF02': 'X', 'KGCF03': 'X', 'KGCF04': 'X', 'KGCF05': 'X', 'KGCF06': 'X', 'KGCF07': 'X', 'KGCF10': 'X', 'KGCF11': 'X', 'KGCF12': 'X', 'KGCF13': 'X', 'KGCF14': 'X', 'KGCF15': 'X', 'KGCF16': 'X', 'KGCF17': 'X', 'KGCF18': 'X','KGCF19': 'X', 'KGCF21': 'X', 'KGCF22': 'X', 'KGCF24': 'X', 'KGCF25': 'X','KGCF26': 'X', 'KGCF27': 'X', 'KGCF28': 'X', 'KGCF29': 'X', 'KGCF30': 'X', 'KGCF32': 'X', 'KGCF33': 'X', 'KGCF34': 'X', 'KGCF35': 'X', 'KGCF36': 'X','KGCF37': 'X', 'KGCF38': 'X', 'KGCF39': 'X', 'KGCF40': 'X', 'KGCF41': 'X','KGCF42': 'X', 'KGCF43': 'X', 'KGCF44': 'X', 'KGCF45': 'X', 'KGCF46': 'X', 'KGCF47': 'X', 'KGCF48': 'X', 'KGCF49': 'X', 'KGCF50': 'X', 'KGCF51': 'X', 'KGCF52': 'X', 'KGCF53': 'X', 'KGCF55': 'X', 'KGCF56': 'X', 'KGCF57': 'X', 'KGCF58': 'X'}

# import file with meta data
metadata_0 = pd.read_csv('original_files/metadata.csv', delimiter=';')
# convert character columns to date objects
metadata_0['birth_date'] = pd.to_datetime(metadata_0['birth_date'])
metadata_0['sample_taken'] = pd.to_datetime(metadata_0['sample_taken'])

for old, new in rename_ids.items():
    metadata_0['identifier'] = metadata_0['identifier'].str.replace(old, new, regex=False)
    
for old, new in rename_patient.items():
    metadata_0['patient'] = metadata_0['patient'].str.replace(old, new, regex=False)

# calculate age in weeks
age_in_weeks = []
for index, row in metadata_0.iterrows():
  age_row = weeks_between(row['birth_date'], row['sample_taken'])
  age_in_weeks.append(age_row)

metadata_0['age_in_weeks'] = age_in_weeks

# keep only samples with patient's age in weeks < 209 (< 4 years)
metadata_1 = metadata_0[metadata_0['age_in_weeks'] < 209] # 2 years  
# store remaining identifiers in list
samples_kept = metadata_1['identifier']

# calculate age in weeks
age_in_weeks = []
for index, row in metadata_0.iterrows():
  age_row = weeks_between(row['birth_date'], row['sample_taken'])
  age_in_weeks.append(age_row)

metadata_0['age_in_weeks'] = age_in_weeks

# keep only samples with patient's age in weeks < 209 (< 4 years)
metadata_1 = metadata_0[metadata_0['age_in_weeks'] < 209] # 2 years  
# store remaining identifiers in list
samples_kept = metadata_1['identifier']

# import file with count data (species)
count_df_0 = pd.read_csv('original_files/count_table_vs_01.csv', delimiter=';')
# keep only species information
count_df_sp_0 = count_df_0.drop(columns=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus'])
# clean column names
new_column_names = count_df_sp_0.columns.str.split('_').str[0]
# add new column names
count_df_sp_1 = count_df_sp_0.set_axis(new_column_names, axis=1, inplace=False)
# remove duplicate rows
count_df_sp_2 = count_df_sp_1.drop_duplicates(subset=None, keep='first')
# set species as row names
count_df_sp_3 = count_df_sp_2.set_index('Species')
# split count data into negative controls and patient samples
count_df_neg_controls = count_df_sp_3.filter(regex=r'(Blank|WA)')
count_df_sp_4 = count_df_sp_3.filter(regex=r'(KGCF|MCF)')
count_df_sp_5 = count_df_sp_4
# store species in list
species_kept = list(count_df_sp_4.index)
samples_kept = list(count_df_sp_4.columns)

# perform clr-transformation on raw count data
count_df_sp_clr_0 = clr(count_df_sp_4.transpose()+1)
count_df_sp_clr_1 = pd.DataFrame(data=count_df_sp_clr_0, index=samples_kept, columns=species_kept)
# re-transpose data frame
count_df_sp_clr_2 = count_df_sp_clr_1.transpose()
# subset data frame
count_df_sp_clr_3 = count_df_sp_clr_2[count_df_sp_clr_2.columns.intersection(samples_kept)]
count_df_sp_clr_4 = count_df_sp_clr_3.rename(columns=rename_ids, inplace=False)
raspir_species = count_df_sp_clr_4.index

# import file with count data (phylum)
# keep only phylum information
count_df_phylum_0 = count_df_0.drop(columns=['Domain', 'Class', 'Order', 'Family', 'Genus', 'Species_merged', "KGCF28_S3"])
# clean column names
new_column_names_phylum = count_df_phylum_0.columns.str.split('_').str[0]
# add new column names
count_df_phylum_1 = count_df_phylum_0.set_axis(new_column_names_phylum, axis=1, inplace=False)
# remove duplicate rows
count_df_phylum_3 = count_df_phylum_1.groupby(['Phylum']).sum()
# split count data into negative controls and patient samples
count_df_neg_controls = count_df_phylum_3.filter(regex=r'(Blank|WA)')
count_df_phylum_4 = count_df_phylum_3.filter(regex=r'(KGCF|MCF)')
count_df_phylum_5 = count_df_phylum_4
# store phyla in list
species_kept_phylum = list(count_df_phylum_4.index)
samples_kept_phylum = list(count_df_phylum_4.columns)
# perform clr-normalisation on raw count data
count_df_phylum_clr_0 = clr(count_df_phylum_4.transpose()+0.01)
count_df_phylum_clr_1 = pd.DataFrame(data=count_df_phylum_clr_0, index=samples_kept_phylum, columns=species_kept_phylum)

# re-transpose data frame
count_df_phylum_clr_2 = count_df_phylum_clr_1.transpose()
# subset data frame
count_df_phylum_clr_3 = count_df_phylum_clr_2[count_df_phylum_clr_2.columns.intersection(samples_kept_phylum)]
count_df_phylum_clr_4 = count_df_phylum_clr_3.rename(columns=rename_ids, inplace=False)

# import file with growth rates
growth_file_00 = pd.read_csv('original_files/growth_rates_vs_03.csv', delimiter=';')
# convert long to wide format
growth_file_01 = growth_file_00.pivot(index='Species', columns='Name', values='Growth_class')
# replace NaNs with string
growth_file_02 = growth_file_01.fillna('unavailable', inplace=False)
growth_file_03 = growth_file_02[growth_file_02.columns.intersection(samples_kept)]
boolean_series = growth_file_03.index.isin(species_kept)
growth_file_04 = growth_file_03[boolean_series]
growth_file_04_1 = growth_file_04[growth_file_04.columns.intersection(samples_kept)]
growth_file_05 = growth_file_04_1.replace('unavailable', 0)
growth_file_06 = growth_file_05.replace('failed', 0)
growth_file_07 = growth_file_06.replace('no growth', 1)
growth_file_08 = growth_file_07.replace('slow', 2)
growth_file_09 = growth_file_08.replace('moderate', 3)
growth_file_10 = growth_file_09.replace('fast', 4)

# add column with row sum
growth_file_10['Total'] = growth_file_10.sum(axis=1)
growth_file_11 = growth_file_10[growth_file_10['Total'] > 10] 
growth_file_12 = growth_file_11.drop(columns='Total')
growth_file_13 = growth_file_12.rename(columns=rename_ids, inplace=False)
growth_file_14 = growth_file_13[growth_file_13.index.isin(raspir_species)]

# import file with growth rates
growth_file_quant_00 = pd.read_csv('original_files/growth_rates_vs_02.csv', delimiter=';')
growth_file_quant_01 = growth_file_quant_00[growth_file_quant_00['Growth_class'] != "failed"]
growth_file_quant_02 = growth_file_quant_01.drop(columns=["No_Reads", "Initial_Bins", "Used_Bins", "Fit_Err"])
growth_file_quant_03 = growth_file_quant_02.set_index('Name')
growth_file_quant_04 = growth_file_quant_03.rename(index=rename_ids, inplace=False)

# export growth file tables
growth_file_13.to_csv("python_generated_files/growth_file_categorical.csv", sep=";", header=True, index=True)
growth_file_quant_04.to_csv("python_generated_files/growth_file_numerical.csv", sep=";", header=True, index=True)

# export metadata
metadata_0.to_csv("python_generated_files/metadata_0.csv", sep=";", header=True, index=False)

# export count data tables
count_df_sp_clr_4.to_csv("python_generated_files/count_df_sp_clr_4.csv", sep=";", header=True, index=True)
count_df_phylum_clr_4.to_csv("python_generated_files/count_df_phylum_clr_4.csv", sep=";", header=True, index=True)
count_df_neg_controls.to_csv("python_generated_files/count_df_neg_controls_clr.csv", sep=";", header=True, index=True)

session_info.show()
