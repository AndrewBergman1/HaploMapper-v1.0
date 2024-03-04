
import pandas as pd
from tabulate import tabulate
import numpy as np

#%% POLISH DATA

# Load the data
df = pd.read_csv('/home/andrewbergman/courses/binp29/pop_gen/01_data/v54.1_1240K_public.anno', sep="\t", low_memory=False)

# Replace commas with dots for the entire DataFrame
df = df.replace(',', '.', regex=True)

# Replace ".." with NaN for proper NA handling
df = df.replace('..', pd.NA, regex=False)

# Convert 'Lat.' to numeric, coercing errors to NaN
df['Lat.'] = pd.to_numeric(df['Lat.'], errors='coerce')
df['Long.'] = pd.to_numeric(df['Long.'], errors='coerce')

# Drop rows where 'Lat.' is NaN
df.dropna(subset=['Lat.'], inplace=True)
df.dropna(subset=['Long.'], inplace=True)
df.dropna(subset=['Y haplogroup (manual curation in ISOGG format)'], inplace=True)

#%% GROUP-BY DATE BINS AND COUNTRY

# Define bins and labels for categorization
bins = [0, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10500, 11000, 11500, 12000]
labels = ['0-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-3500', '3501-4000', '4001-4500', '4501-5000', '5001-5500', '5501-6000', '6001-6500', '6501-7000', '7001-7500', '7501-8000', '8001-8500', '8501-9000', '9001-9500', '9501-10000', '10001-10500', '10501-11000', '11001-11500', '11501-12000']



# Categorize observations into bins
df['DateBins'] = pd.cut(df['Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'], bins=bins, labels=labels, right=False)

# Convert 'Locality' and 'DateBins' to categorical
#df['Locality'] = df['Locality'].astype('category')
#df['DateBins'] = df['DateBins'].astype('category')

df['CombinedBins'] = df['Locality'] + " ( " + df['DateBins'].astype('str') + " ) "
#df['CombinedBins'] = df['CombinedBins'].astype('category')

#%% Y HAPLOTYPE DUMMY VARIABLES

# Y HAPLOTYPE DUMMY VARIABLES
y_dummies = pd.get_dummies(df['Y haplogroup (manual curation in ISOGG format)'])
cols_to_keep = [col for col in y_dummies.columns if not col.startswith('n/a') and col[0].isalpha()]
y_dummies = y_dummies[cols_to_keep].copy()

unique_starting_letters = set(col[0] for col in y_dummies.columns)
for letter in unique_starting_letters:
    cols_to_sum = [col for col in y_dummies.columns if col.startswith(letter)]
    y_dummies[f'{letter}_sum'] = y_dummies[cols_to_sum].sum(axis=1)

y_dummies_summed = y_dummies.filter(like='_sum').copy()
y_dummies_summed["y_haplos_sum"] = y_dummies_summed.sum(axis=1)

# Correct Division Process
# Convert all columns to float for division
y_dummies_summed = y_dummies_summed.astype(float)

# Append the modified dummies to the original DataFrame and proceed with aggregation
df = pd.concat([df, y_dummies_summed], axis=1)


#%% CREATE TABLE
# Sort by 'Locality' and 'DateBins'
df = df.sort_values(by=['CombinedBins'])

# Define the aggregation dictionary, initially for 'Lat.' and 'Long.'
agg_dict = {'Lat.': 'mean', 'Long.': 'mean'}

# Extend the aggregation dictionary to include sums for all Y haplogroup dummy variables
for col in y_dummies_summed.columns:
    if col.endswith('_sum'):
        agg_dict[col] = 'sum'

# Aggregate to find mean 'Lat.' for each group
aggregated_df = df.groupby(['CombinedBins'], as_index=False, observed=True).agg(agg_dict)


# Print the first row of the aggregated DataFrame
#print(tabulate(aggregated_df.head(20), headers='keys', tablefmt='psql'))

#print("Number of columns:", aggregated_df.shape[1])
#print(aggregated_df["y_haplos_sum"])

#print(aggregated_df)

# Ensure 'y_haplos_sum' is numeric
aggregated_df['y_haplos_sum'] = pd.to_numeric(aggregated_df['y_haplos_sum'], errors='coerce')

# Perform division only on numeric columns
numeric_cols = aggregated_df.select_dtypes(include=[np.number]).columns.drop('y_haplos_sum')

for col in numeric_cols:
    if col != "Lat." and col != "Long." :
        aggregated_df[col] = aggregated_df[col].div(aggregated_df['y_haplos_sum'])*100

# Handle any potential NaN values resulting from the division
aggregated_df.fillna(0, inplace=True)

with open('/home/andrewbergman/courses/binp29/pop_gen/table.txt', 'w') as f:
    f.write(tabulate(aggregated_df, headers='keys', tablefmt='psql'))
