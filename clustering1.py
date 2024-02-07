import pandas as pd
import warnings
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import plotly.express as px
import random

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

file_path = 'C:\\compsci\\4150\\clustering-1\\GSE64881_segmentation_at_30000bp.passqc.multibam.txt'

# Read the text file into a DataFrame
df = pd.read_csv(file_path, delimiter='\t')

filtered_data = df[(df['chrom'] == 'chr13') & (df['start'] > 21700000) & (df['stop'] < 24100000)]
# Ignore 0 only columns & first 3 columns
ignored_columns = filtered_data.iloc[:, 3:].loc[:, (filtered_data.iloc[:, 3:] != 0).any()]

num_rows, num_columns = filtered_data.shape

total_ones_col = ignored_columns.sum()

# Average windows present in an NP
avg_windows = total_ones_col.mean()

# Least and Most windows present in an NP
least_windows_id = total_ones_col.idxmin()
least_windows_count = total_ones_col.min()
most_windows_id = total_ones_col.idxmax()
most_windows_count = total_ones_col.max()
# Count the number of rows that contain at least one "1"
num_rows_with_ones = ignored_columns.apply(lambda row: 1 if 1 in row.values else 0, axis=1).sum()
# Calculate the average number of rows with at least one "1"
average_rows_with_ones = num_rows_with_ones / num_rows
num_nps = 0
for column in ignored_columns:
    column_values = filtered_data[column].sum()
    if column_values > 0:
        num_nps += 1
 
print(f"Num windows: {num_rows}")
print(f"Num NPs: {num_nps}")
print(f"Average Windows Present in an NP: {avg_windows}")
print(f"NP: {least_windows_id} has the smallest amount of windows with {least_windows_count} window(s)")
print(f"NP: {most_windows_id} has the largest amount of windows with {most_windows_count} window(s)")
print(f"Average number of NPs with a window detected in H1: {average_rows_with_ones}")

# Columns to ignore
columns_to_ignore = df.columns[:3]

# get rid of empty columns and sum 1s
non_empty_columns = df.drop(columns=columns_to_ignore).columns
sum_of_ones = df[non_empty_columns].sum()

# Calculate bin edges and get thresholds for columns
thresholds = pd.cut(sum_of_ones, bins=5, retbins=True)[1]

# Calculate bin edges and get thresholds for rows
row_thresholds = pd.cut(filtered_data.iloc[:, 3:].sum(axis=1), bins=10, retbins=True)[1]

# Define bins for scale mapping
bins = [-float('inf')] + thresholds.tolist() + [float('inf')]
row_bins = [-float('inf')] + row_thresholds.tolist() + [float('inf')]

# DataFrame to store scale values for each column
scale_df = pd.DataFrame()
row_scale_df = pd.DataFrame()
comp_pos = 0
for index, row in filtered_data.iterrows():
    row_values = row.iloc[3:]  # Assuming the first 3 columns are not part of the scaling

    # Skip the row if it is empty or has sum of 0
    if row_values.empty or row_values.sum() == 0:
        continue

    sum_of_ones_row = row_values.sum()

    # Assign scale value to the row 
    row_scale_value = pd.cut([sum_of_ones_row], bins=row_bins, labels=False, right=False).astype(int)[0]
    # To total sum comp position to eventually find average
    comp_pos = comp_pos + row_scale_value

    # Assign the scale value to the row
    row_scale_df.loc[index, 'row_scale'] = row_scale_value

rad_pos = 0
# Iterate through columns in the filtered_data
for column in filtered_data.columns[3:]:
    column_values = filtered_data[column].dropna()

    # Skip the column if empty or has sum of 0
    if column_values.empty or column_values.sum() == 0:
        continue

    sum_of_ones_col = column_values.sum()
    

    # Use cut to categorize into bins and assign scale values
    scale_value = pd.cut([sum_of_ones_col], bins=bins, labels=False, right=False).astype(int) 
    rad_pos = rad_pos + scale_value
    #print(f'{column} - Times Found: {sum_of_ones_col}, Scale Value: {scale_value[0]}')

    # Assign the scale value to each row in the column
    scale_df[column + '_scale'] = scale_value[0]

print(f"Average radial position for H1 is: ", (rad_pos / num_columns))
print(f"Average compaction  for H1 is: ", (comp_pos / num_rows))

num_columns = len(ignored_columns.columns)
jaccard_matrix = np.zeros((num_columns, num_columns))

for i in range(num_columns - 1):
    column_i = ignored_columns.iloc[:, i]
    for j in range(i + 1, num_columns):
        column_j = ignored_columns.iloc[:, j]

        intersection = sum((column_i == 1) & (column_j == 1))
        #New Denominator
        
        usumi = sum((column_i == 1))
        usumj = sum((column_j == 1))
        union = min(usumi, usumj)

        jaccard_index = intersection / union

        # upper triangular part of the matrix
        jaccard_matrix[i, j] = jaccard_index

# Copy 
jaccard_matrix = jaccard_matrix + jaccard_matrix.T - np.diag(jaccard_matrix.diagonal())

# Calculate Jaccard distance matrix
jaccard_distance_matrix = 1 - jaccard_matrix
random_num = random.sample(range(1,162),3)
print(random_num)
#Create cluster lists
c1 = []
c2 = []
c3 = []

for i in range(1,162):
    c1dist = jaccard_distance_matrix[random_num[0],i]
    c2dist = jaccard_distance_matrix[random_num[1],i]
    c3dist = jaccard_distance_matrix[random_num[2],i]
    assignment = min(c1dist,c2dist,c3dist)
    if assignment == c1dist:
        c1.append(i)
    if assignment == c2dist:
        c2.append(i)
    if assignment == c3dist:
        c3.append(i)
    
print(f"Cluster {random_num[0]}: ",c1)    
print(f"Cluster {random_num[1]}: ",c2)
print(f"Cluster {random_num[2]}: ",c3)
        
    
        
    
    

# Create DataFrames for Jaccard index and distance
jaccard_df = pd.DataFrame(jaccard_matrix, columns=ignored_columns.columns, index=ignored_columns.columns)
jaccard_distance_df = pd.DataFrame(jaccard_distance_matrix, columns=ignored_columns.columns, index=ignored_columns.columns)

# Save the DataFrames to CSV files
jaccard_df.to_csv('jaccard_matrix.csv')
jaccard_distance_df.to_csv('jaccard_distance_matrix.csv')

# Specify the file paths in the current directory
fpath_jaccard = os.path.join('C:\\compsci\\4150\\clustering-1', 'jaccard_matrix.csv')
fpath_distance = os.path.join('C:\\compsci\\4150\\clustering-1', 'jaccard_distance_matrix.csv')

# Save the DataFrames to CSV files in the current directory
jaccard_df.to_csv(fpath_jaccard)
jaccard_distance_df.to_csv(fpath_distance)

# Print or use the matrices as needed
print("Jaccard Matrix:")
print(jaccard_matrix)
print("\nJaccard Distance Matrix:")
print(jaccard_distance_matrix)
    
fig = px.imshow(jaccard_matrix, text_auto=True)
fig.write_html('jaccard-heat.html')

fig2 = px.imshow(jaccard_distance_matrix, text_auto=True)
fig2.write_html('jaccard-heat1.html')


# Plot the scatter plot
#plt.scatter(range(len(sum_of_ones)), sum_of_ones)
#plt.title('Times NP was Found in a Slice')
#plt.xlabel('Slice Index')
#plt.ylabel('Times NP Was Found')
#plt.show()
