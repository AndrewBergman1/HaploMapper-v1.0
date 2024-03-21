'''
Run the program as follows: Python HaploMapper.py

You will be prompted with the paths to the files necessary for the program to run. All of the download links are present in the provided ReadMe file.

'''

import pandas as pd
from tabulate import tabulate
import numpy as np
import plotly.express as px
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import re
import os
from flask_caching import Cache

def boot():

    '''
    Prompts the user for file paths required to run HaploMapper.

    output: File paths.

    '''

    print("Welcome to HaploMapper! Before you continue, ensure that you have downloaded the files from the ReadMe.")
    
    mt_snp_file = input("Provide the path to mitochondrial SNP file (e.g., /path/to/mtSNP.csv): ")
    mt_basal_file = input("Provide the path to the mitochondrial haplogroup phylogeny file (e.g., /path/to/mtPhylo.csv): ")
    y_snp_file = input("Provide the path to yDNA SNP file (e.g., /path/to/ySNP.csv): ")
    y_locus_file = input("Provide the path to the locus file (e.g., /path/to/locusFile.csv): ")
    y_basal_file = input("Provide the path to the yDNA haplogroup phylogeny file (e.g., /path/to/yPhylo.csv): ")
    anno_file = input("Provide the path to the annotation file (e.g., /path/to/anno.csv): ")
    bin_choice = input("Enter the time resolution of interest (e.g., 100 -> 100 year resolution on map, 1000 -> 1000 year resolution on map): ")

    # Optionally, validate that each provided path exists and is a file
    for file_path in (mt_snp_file, mt_basal_file, y_snp_file, y_locus_file, y_basal_file, anno_file):
        if not os.path.isfile(file_path):
            print(f"Error: The provided path does not exist or is not a file: {file_path}")
            return "", "", "", "", "", "", int(bin_choice) #CHANGE TO NONE WHEN HANDING IN!!! 

    return mt_snp_file, mt_basal_file, y_snp_file, y_locus_file, y_basal_file, anno_file, int(bin_choice)

def open_data(anno_file) :

    '''
    Reads data from the provide file paths.

    input: AADR dataset
    output: dataframes for Y-chromosome and mitochondrial DNA.
    '''
    # Load the data
    df = pd.read_csv(anno_file, sep="\t", low_memory=False)

    # Replace commas with dots for the entire Dataframe
    df = df.replace(',', '.', regex=True)

    # Replace ".." with NaN 
    df = df.replace('..', pd.NA, regex=False)

    df.dropna(subset=['Lat.'], inplace=True)
    df.dropna(subset=['Long.'], inplace=True)

    # Remove non-numeric sample ages
    sample_ages = [str(age) if str(age).isnumeric() else 'n/a' for age in df['Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]']]
    df['Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'] = pd.to_numeric(sample_ages, errors='coerce')
    
    # Convert 'Lat.' to numeric, coercing errors to NaN
    df['Lat.'] = pd.to_numeric(df['Lat.'], errors='coerce')
    df['Long.'] = pd.to_numeric(df['Long.'], errors='coerce')

    y_df = df
    mt_df = df

    # Drop rows where haplogroups are NaN for both dataframes
    y_df.dropna(subset=['Y haplogroup (manual curation in ISOGG format)'], inplace=True)
    mt_df.dropna(subset=['mtDNA haplogroup if >2x or published'], inplace=True)
    # Return dataframes
    return y_df, mt_df

def create_bins(y_df, mt_df, bin_choice) :
    '''
    Bin the data frames based on political entity (country) and sample age.

    input: dataframes for Y-chromosome and mitochondrial DNA.
    output: binned dataframes for Y-chromosome and mitochondrial DNA.
    '''
    # Bin size is defined to be compatible with future additions to AADR dataset
    start = 0
    end = 120000
    bin_size = bin_choice # Bin choice is defined by the user.

    # Create bins and labels based on the start and end variables.
    bins = list(range(start, end + bin_size, bin_size))
    labels = [f"{i}-{i + bin_size - 1}" for i in range(start, end, bin_size)]

    # Categorize observations into bins    
    mt_df['DateBins'] = pd.cut(mt_df['Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'], bins=bins, labels=labels, right=False)
    y_df['DateBins'] = pd.cut(y_df['Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'], bins=bins, labels=labels, right=False)

    y_df['CombinedBins'] = y_df['Political Entity'] + " (" + y_df['DateBins'].astype('str') + "BP)"
    mt_df['CombinedBins'] = mt_df['Political Entity'] + " (" + mt_df['DateBins'].astype('str') + "BP)"

    # Return dataframes
    return y_df, mt_df

def findYBasalHaplogroups(y_df, y_snp_file, y_basal_file, y_locus_file):

    '''
    Identifies the basal haplogroups from the haplogroups provided in the AADR dataset.

    input: binned dataframes for Y-chromosome and mitochondrial DNA as well as Y_SNP_FILE, Y_LOCUS_FILE and Y_PHYLO_FILE
    output: Y-chromosome dataframe, now containing the basal haplogroups.
    '''
    # regular expression for replacing spaces with tabs
    space_to_tab = re.compile(r'\s+')

    # Function that loads a dictionary from a file
    # key_index and value_index correspond to the fields where they are found, respectively.
    def read_dict_from_file(file_path, key_index, value_index, split_char='\t'):
        result_dict = {}
        with open(file_path, 'r') as file:
            for line in file:
                parts = space_to_tab.sub(split_char, line.strip()).split(split_char)
                if len(parts) >= max(key_index, value_index) + 1:
                    result_dict[parts[key_index].strip()] = parts[value_index].strip()
        return result_dict

    # Call the read_dict_from_file function.
    y_basal_dict = read_dict_from_file(y_basal_file, 0, 1)
    y_snp_dict = read_dict_from_file(y_snp_file, 1, 3)
    y_locus_dict = read_dict_from_file(y_locus_file, 1, 3)

    # Clean the data using list comprehension
    y_haplo_data = [haplo.split(";")[0].split("-")[0].split("(")[0].replace('~', "").replace('*', "") for haplo in y_df["Y haplogroup (manual curation in ISOGG format)"].tolist()]

    # Reverse the dictionary for easier use in the loop below
    y_snp_dict_inv = {v: k for k, v in y_snp_dict.items()}

    # Update haplogroup data based on dictionaries
    loop = 0

    # Loop as long as there are haplogroups that aren't references by one letter (that isnt n/a(Female)), Max 30 loops
    while any(len(item) > 1 and item != 'n/a(Female)' for item in y_haplo_data) and loop < 50:
        loop += 1
        # Loop through every haplogroup
        for index, haplo in enumerate(y_haplo_data):
            # If the haplogroup is refered to as with more than a character 
            if len(haplo) > 1:
                # Check phylo file
                if haplo in y_basal_dict and y_basal_dict[haplo] not in ['#', "Root"]:
                    y_haplo_data[index] = y_basal_dict[haplo]
                # Check SNP file
                elif haplo in y_snp_dict_inv:
                    y_haplo_data[index] = y_snp_dict_inv[haplo]
                # Check locus file
                elif haplo in y_locus_dict:
                    y_haplo_data[index] = y_locus_dict[haplo]
                # Replace haplogroup with n/a if its not found
                else:
                    y_haplo_data[index] = "n/a"
        # Replace the column in the y-chromosome dataframe w/ y_haplo_data
        y_df["Y haplogroup (manual curation in ISOGG format)"] = y_haplo_data
    
    # Return dataframe, now containing the basal haplogroups
    return y_df

def findMTBasalHaplogroups(mt_df, mt_snp_file, mt_basal_file):
    '''
    Identifies the basal haplogroups from the haplogroups provided in the AADR dataset.

    input: binned dataframes for Y-chromosome and mitochondrial DNA as well as Y_SNP_FILE, Y_LOCUS_FILE and Y_PHYLO_FILE
    output: mitochondrial DNA dataframe, now containing the basal haplogroups.
    '''
    # Function that loads a dictionary from a file
    # key_index and value_index correspond to the fields where they are found, respectively.
    def read_dict_from_file(file_path, key_index, value_index, split_char='\t'):
        result_dict = {}
        with open(file_path, 'r') as file:
            for line in file:
                parts = space_to_tab.sub(split_char, line.strip()).split(split_char)
                if len(parts) >= max(key_index, value_index) + 1:
                    result_dict[parts[key_index].strip()] = parts[value_index].strip()
        return result_dict
    
    # regular expression for replacing spaces with tabs
    space_to_tab = re.compile(r'\s+')

    # Create mtDNA haplogroup dictionary from a phylogenetic tree file
    mt_basal_dict = read_dict_from_file(mt_basal_file, 0, 1)

    # Create mutant dictionary from SNP file (keeping the provided file path even though it's labeled "y")
    mutant_dict = read_dict_from_file(mt_snp_file, 1, 3)

    # Clean the data using list comprehension
    mt_haplo_data = mt_df['mtDNA haplogroup if >2x or published'].tolist()
    mt_haplo_data = [haplo.split(";")[0].split("-")[0].split("(")[0].replace('~', "").replace('*', "") for haplo in mt_haplo_data]

    # Reverse the dictionary for easier use in the loop below
    mutant_dict_inv = {v: k for k, v in mutant_dict.items()}

    # Update haplogroup data based on dictionaries
    loop = 0
    # Loop as long as there are haplogroups that aren't references by one letter (that isnt n/a(Female)), Max 30 loops
    while any(len(item) > 1 and item != 'n/a(Female)' for item in mt_haplo_data) and loop < 50:
        loop += 1

        # Loop through every haplogroup
        for index, haplo in enumerate(mt_haplo_data):
            # If the haplogroup is refered to as with more than a character 
            if len(haplo) > 1:  
                # check phylo file
                if haplo in mt_basal_dict and mt_basal_dict[haplo] not in ['#', "Root"]:
                    mt_haplo_data[index] = mt_basal_dict[haplo]
                # check mutant file
                elif haplo in mutant_dict_inv:
                    mt_haplo_data[index] = mutant_dict_inv[haplo]
                # replace with n/a if not present
                else:
                    mt_haplo_data[index] = "n/a"
    
    # Update data frame with basal haplogroups
    mt_df['Updated mtDNA haplogroup'] = mt_haplo_data  

    # Return dataframe, now containing basal haplogroups
    return mt_df

def createDummyVariables (y_df, mt_df) :
    '''
    Creates dummy variables for each haplogroup. Each dummy varaible carries a binary representation of whether it occurs in a sample or not. 
    
    input: mitochondrial DNA dataframe and Y-chromosome dataframe containing basal haplogroups.
    output: mitochondrial DNA dataframe and Y-chromosome dataframe containing basal haplogroups and dummy-variable dataframe for each.
    '''
    # Y HAPLOTYPE DUMMY VARIABLES
    y_dummies = pd.get_dummies(y_df['Y haplogroup (manual curation in ISOGG format)']) # The y haplogroups
    mt_dummies = pd.get_dummies(mt_df['mtDNA haplogroup if >2x or published']) # The mt haplogroups

    # Keep columns that dont start with 'n', as they are n/a or nan
    y_cols_to_keep = [col for col in y_dummies.columns if not col.startswith('n') and col[0].isalpha()]
    mt_cols_to_keep = [col for col in mt_dummies.columns if not col.startswith('n') and col[0].isalpha()]

    # Subset the dummy variables dataframes
    y_dummies = y_dummies[y_cols_to_keep].copy()
    mt_dummies = mt_dummies[mt_cols_to_keep].copy()

    # Some haplogroups cannot be further tracked than 'I2', to which the basal haplogroup is I. Since HaploMapper visualises the basal haplogroups, the first letter is retrieved
    y_unique_starting_letters = set(col[0] for col in y_dummies.columns)
    mt_unique_starting_letters = set(col[0] for col in mt_dummies.columns)

    # Identify the haplogroup-carrying columns
    for letter in y_unique_starting_letters:
        cols_to_sum = [col for col in y_dummies.columns if col.startswith(letter)]
        y_dummies[f'{letter}_y_sum'] = y_dummies[cols_to_sum].sum(axis=1)
   
    # Identify the haplogroup-carrying columns
    for letter in mt_unique_starting_letters:
        cols_to_sum = [col for col in mt_dummies.columns if col.startswith(letter)]
        mt_dummies[f'{letter}_mt_sum'] = mt_dummies[cols_to_sum].sum(axis=1)

    # Copy the haplogroup sum. This corresponds to the number of individuals found at any given sampling site 
    y_dummies_summed = y_dummies.filter(like='_sum').copy()
    y_dummies_summed["y_haplos_sum"] = y_dummies_summed.sum(axis=1)

    # Copy the haplogroup sum. This corresponds to the number of individuals found at any given sampling site 
    mt_dummies_summed = mt_dummies.filter(like='_sum').copy()
    mt_dummies_summed["mt_haplos_sum"] = mt_dummies_summed.sum(axis=1)

    # Convert all columns to float for division
    y_dummies_summed = y_dummies_summed.astype(float)
    mt_dummies_summed = mt_dummies_summed.astype(float)

    # Append the modified dummies to the original DataFrame and proceed with aggregation
    y_df = pd.concat([y_df, y_dummies_summed], axis=1)
    mt_df = pd.concat([mt_df, mt_dummies_summed], axis=1)

    return y_df, mt_df

def createTable(combined_df) :
    '''
    Concatenates the dummy table and the dataframes correpsodning to Y-chromosome and mitochondrial DNA AADR.
    
    input: mitochondrial DNA dataframe and Y-chromosome dataframe containing basal haplogroups and dummy-variable dataframe for each.
    output: A data frame per genetic carrier which is composed of the dummy variables dataframe concatenated with the corresponding longitudes, latitudes and bins.
    '''
    def custom_aggregate(series):
        if series.dtype == 'float64' or series.dtype == 'int64':
            return series.sum()  # Sum for numeric columns
        else:
            # Return the first non-NA string value in the series
            return series.dropna().iloc[0] if not series.dropna().empty else None


    # Group by long, lat and combinedbins in order to visualise overlapping data points 
    combined_df_inds = combined_df.groupby(['Lat.', 'Long.', 'CombinedBins']).agg(custom_aggregate).reset_index()

    # Define the aggregation dictionary
    agg_dict = {}
    
    # Include all columns carrying haplogroup information
    for col in combined_df.columns:
        if col.endswith('_sum') or col.endswith('y_haplos_sum') or col.endswith('mt_haplos_sum'):
            agg_dict[col] = 'sum'
    
    # Aggregate the observations to the nation-time bins (ex Sweden 0-999 BP)
    combined_df = combined_df.groupby(['CombinedBins'], as_index=False, observed=True).agg(agg_dict)
    
    # Ensure the sum of haplogroups are numeric
    combined_df['y_haplos_sum'] = pd.to_numeric(combined_df['y_haplos_sum'], errors='coerce')
    combined_df['mt_haplos_sum'] = pd.to_numeric(combined_df['mt_haplos_sum'], errors='coerce')

    # Return a data frame, binned by nation and time. For each bin, the y-haplogroups and mt-haplogroups are registered 
    return combined_df, combined_df_inds

def saveToFile(combined_df_inds, combined_df_nations) :
    '''
    Saves the frequency tables.

    input: A data frame per genetic carrier which is composed of the dummy variables dataframe concatenated with the corresponding longitudes, latitudes and bins.
    '''
    # Save the combined_df_inds as a CSV file
    combined_df_inds.to_csv('/home/andrewbergman/courses/binp29/pop_gen/all_observations.csv', index=False)

    # Save the combined_df_nations as a CSV file
    combined_df_nations.to_csv('/home/andrewbergman/courses/binp29/pop_gen/national_obersvations.csv', index=False)

def createWebApplication(combined_df_inds, combined_df_nations) :
    '''
    Creates a dash application. Interactive map is generated by px.scatterplot() and pie charts are generated by px.pie(). display_y_click_data and display_mt_click_data provide w/ corr. callback functions provide the interactivity with the map. 

    input: Dataframes containing basal haplogroups, latitude, longitude and bin.
    output: Dash application.
    '''
    app = dash.Dash(__name__)
    
    # Color map used to standardise the haplogroup colours
    color_map = {
        'A': '#3182bd',
        'B': '#3182bd',
        'C': '#6baed6',
        'D': '#9ecae1',
        'E': '#c6dbef',
        'F': '#e6550d',
        'G': '#e6550d',
        'H': '#fd8d3c',
        'I': '#fdae6b',
        'J': '#fdd0a2',
        'K': '#31a354',
        'L': '#31a354',
        'M': '#74c476',
        'N': '#a1d99b',
        'O': '#c7e9c0',
        'P': '#756bb1',
        'Q': '#756bb1',
        'R': '#9e9ac8',
        'S': '#bcbddc',
        'T': '#dadaeb',
        'U': '#636363',
        'V': '#636363',
        'W': '#969696',
        'X': '#bdbdbd',
        'Y': '#d9d9d9',
        'Z': '#d9d9d9'
    }

    # Bin with respect to time
    timeBins = combined_df_inds['CombinedBins'].str.split("(", n=1, expand = True)[1].str.split(")", n=1, expand=True)[0] 
    parts = combined_df_inds['CombinedBins'].str.split(" ", n=1, expand = True)
    
    # Retrieve the countries (political entities) for color-coding the interactive map.
    if len(parts) > 2 : 
        country = combined_df_inds['CombinedBins'].str.split(" ", n=1, expand = True)[0] + " " + combined_df_inds['CombinedBins'].str.split(" ", n=1, expand = True)[1]
    else : 
        country = combined_df_inds['CombinedBins'].str.split(" ", n=1, expand = True)[0]

    # add times and ages to the dataframe
    combined_df_inds['Time'] = timeBins
    combined_df_inds['Country'] = country
    combined_df_inds['SortKey'] = combined_df_inds['Time'].apply(lambda x: int(x.split('-')[0]))

    # Step 2: Sort the dataframe based on the new 'SortKey' column
    combined_df_inds_sorted = combined_df_inds.sort_values(by='SortKey')

    # Now generate the interactive map using the sorted dataframe
    fig = px.scatter_geo(combined_df_inds_sorted, lat='Lat.', lon='Long.',
                        projection="natural earth",
                        custom_data=["Long.", "Lat.", "CombinedBins", 'Country'],
                        color="Time",
                        hover_data=['Country'],
                        template="seaborn")
    
    # Update the layout to show country borders and names
    fig.update_geos(
        showcountries=True, 
        countrycolor="RebeccaPurple", # Color for the country borders
        showcoastlines=True, 
        coastlinecolor="RebeccaPurple", # Color for the coastline
        showland=True, 
        landcolor='LightGreen', # Color for land
        showocean=True, 
        oceancolor='LightBlue', # Color for ocean
        showsubunits=True, 
        subunitcolor="RebeccaPurple", # Color for subunits
        showframe=False 
    )
    # Alter the size of the map in the application
    fig.update_layout(height=550, width=1500)

    # your existing layout
    app.layout = html.Div([
        html.H1("HaploMapper", style={'font-family': 'Arial, sans-serif'}),  # Main title
        html.H2("Welcome to HaploMapper, a tool visualizing Allen Ancient DNA Resource (AADR). The map below hosts map markers, each corresponding to observations of basal haplogroups. The map-markers are color-coded by the samples' age. Upon interacting with a map marker, the distribution of basal haplogroups on the Y-chromosome and the mitochondrial DNA are illustrated in pie charts. The data corresponding to the visualisations are saved on your local machine, in the same directory that you ran HaploMapper from.", style={'font-family': 'Arial, sans-serif'}),  # Description
        dcc.Graph(id='map-plot', figure=fig),  # The map plot
        dcc.Store(id='store-data', data=combined_df_inds.to_dict('records')),

        # Container for pie charts
        html.Div([
            html.Div([
                html.H2("Y Basal Haplogroup Distribution", style={'font-family' : 'Arial, sans-serif'}),  # Description for the Y haplogroup pie chart
                dcc.Graph(id='y_pie_chart'),  # The Y haplogroup pie chart
            ], style={'display': 'inline-block', 'width': '50%'}),
        
            # Container for the MT haplogroup pie chart
            html.Div([
                html.H2("MT Basal Haplogroup Distribution", style={'font-family': 'Arial, sans-serif'}),  # Description for the MT haplogroup pie chart
                dcc.Graph(id='mt_pie_chart'),  # The MT haplogroup pie chart
            ], style={'display': 'inline-block', 'width': '50%'}),
            html.Div([
                html.H2("Y Basal Haplogroup NATION Distribution", style={'font-family': 'Arial, sans-serif'}),  # Description for the Y haplogroup pie chart
                dcc.Graph(id='y_pie_chart_nation'),  # The Y haplogroup pie chart
            ], style={'display': 'inline-block', 'width': '50%'}),
            # Container for the MT haplogroup pie chart
            html.Div([
                html.H2("MT Basal Haplogroup NATION Distribution", style={'font-family': 'Arial, sans-serif'}),  # Description for the MT haplogroup pie chart
                dcc.Graph(id='mt_pie_chart_nation'),  # The MT haplogroup pie chart
            ], style={'display': 'inline-block', 'width': '50%'}),
        ]),
    ],
    style={'textAlign': 'center', 'width': '80%', 'margin': 'auto'})  # Styling

    # Callback function to display y basal haplogroups upon clicking on a map marker
    @app.callback(
        Output('y_pie_chart', 'figure'),
        [Input('map-plot', 'clickData')]
    )
    def display_y_click_data(clickData):
        if clickData:
            longitude = clickData['points'][0]['customdata'][0]
            latitude = clickData['points'][0]['customdata'][1]
            bin = clickData['points'][0]['customdata'][2]

            y_row = combined_df_inds[(combined_df_inds['Long.'] == longitude) & 
                                    (combined_df_inds['Lat.'] == latitude) & 
                                    (combined_df_inds['CombinedBins'] == bin)]

            if not y_row.empty:
                y_row = y_row.iloc[0]
                # Create a DataFrame for haplogroups and their frequencies
                y_haplogroups = pd.DataFrame({
                    'Haplogroup': [col.replace('_y_sum', '') for col in y_row.index if col.endswith('_y_sum') and y_row[col] > 0],
                    'Frequency': [y_row[col] for col in y_row.index if col.endswith('_y_sum') and y_row[col] > 0]
                })

            if not y_haplogroups.empty:
                title_text = f'Basal Y Haplogroup Distribution at the Selected Sampling Site: {longitude}; {latitude}.<br>There are {int(y_row["y_haplos_sum"])} samples at this site.'
                fig = px.pie(y_haplogroups, names='Haplogroup', values='Frequency', title=title_text)
                # Apply the color mapping
                fig.update_traces(marker=dict(colors=[color_map[haplo] for haplo in y_haplogroups['Haplogroup']]))
                return fig

        return px.pie(title='No Y Haplogroup Data Selected')
        
    # Callback function to display mt basal haplogroups upon clicking on a map marker
    @app.callback(
        Output('mt_pie_chart', 'figure'),
        [Input('map-plot', 'clickData')]
    )
    # Function retrieving data from the mt dataframe and returns a pie chart
    def display_mt_click_data(clickData):
        if clickData:
            longitude = clickData['points'][0]['customdata'][0]
            latitude = clickData['points'][0]['customdata'][1]
            bin = clickData['points'][0]['customdata'][2]
            
            mt_row = combined_df_inds[(combined_df_inds['Long.'] == longitude) & 
                                    (combined_df_inds['Lat.'] == latitude) & 
                                    (combined_df_inds['CombinedBins'] == bin)]
            
            if not mt_row.empty:
                mt_row = mt_row.iloc[0]
                mt_haplogroups = pd.DataFrame({
                    'Haplogroup': [col.replace('_mt_sum', '') for col in mt_row.index if col.endswith('_mt_sum') and mt_row[col] > 0],
                    'Frequency': [mt_row[col] for col in mt_row.index if col.endswith('_mt_sum') and mt_row[col] > 0]
                })

                if not mt_haplogroups.empty:
                    title_text = f'Basal MT Haplogroup Distribution at the Selected Sampling Site: {longitude}; {latitude}.<br>There are {int(mt_row["mt_haplos_sum"])} samples at this site.'
                    fig = px.pie(mt_haplogroups, names='Haplogroup', values='Frequency', title=title_text)
                    fig.update_traces(marker=dict(colors=[color_map[haplo] for haplo in mt_haplogroups['Haplogroup']]))
                    return fig

        return px.pie(title='No MT Haplogroup Data Selected')
    
    # Callback function to display y basal haplogroups upon clicking on a map marker
    @app.callback(
        Output('y_pie_chart_nation', 'figure'),
        [Input('map-plot', 'clickData')]
    )
    # Function retrieving data from the Y dataframe and returns a pie chart
    def display_y_nation_click_data(clickData):
        if clickData:
            bin = clickData['points'][0]['customdata'][2]
            y_row_nation = combined_df_nations[combined_df_nations['CombinedBins'] == bin]
            
            if not y_row_nation.empty:
                y_row_nation = y_row_nation.iloc[0]
                y_haplogroups_nation = pd.DataFrame({
                    'Haplogroup': [col.replace('_y_sum', '') for col in y_row_nation.index if col.endswith('_y_sum') and y_row_nation[col] > 0],
                    'Frequency': [y_row_nation[col] for col in y_row_nation.index if col.endswith('_y_sum') and y_row_nation[col] > 0]
                })

                if not y_haplogroups_nation.empty:
                    fig = px.pie(y_haplogroups_nation, names='Haplogroup', values='Frequency',
                                title='Basal Y Haplogroups Distribution in ' + bin + ' There are ' + str(int(y_row_nation['y_haplos_sum'])) + " samples in total.")
                    fig.update_traces(marker=dict(colors=[color_map[haplo] for haplo in y_haplogroups_nation['Haplogroup']]))
                    return fig

        return px.pie(title='No Y Haplogroup Data Selected')
    
    # Callback function to display mt basal haplogroups upon clicking on a map marker
    @app.callback(
        Output('mt_pie_chart_nation', 'figure'),
        [Input('map-plot', 'clickData')]
    )
    # Function retrieving data from the mt dataframe and returns a pie chart
    def display_mt_nation_click_data(clickData):
        if clickData:
            bin = clickData['points'][0]['customdata'][2]

            mt_row_nation = combined_df_nations[combined_df_nations['CombinedBins'] == bin]
            if not mt_row_nation.empty:
                mt_row_nation = mt_row_nation.iloc[0]
                mt_haplogroups_nation = pd.DataFrame({
                    'Haplogroup': [col.replace('_mt_sum', '') for col in mt_row_nation.index if col.endswith('_mt_sum') and mt_row_nation[col] > 0],
                    'Frequency': [mt_row_nation[col] for col in mt_row_nation.index if col.endswith('_mt_sum') and mt_row_nation[col] > 0]
                })

                if not mt_haplogroups_nation.empty:
                    fig = px.pie(mt_haplogroups_nation, names='Haplogroup', values='Frequency',
                                title='Basal MT Haplogroups Distribution in ' + bin + ' There are ' + str(int(mt_row_nation['mt_haplos_sum'])) + " samples in total.")
                    fig.update_traces(marker=dict(colors=[color_map[haplo] for haplo in mt_haplogroups_nation['Haplogroup']]))
                    return fig

        return px.pie(title='No MT Haplogroup Data Selected')
    
    if __name__ == '__main__':
        app.run_server(debug=False, port=8052) 

#Retrive the file paths
mt_mut_file, mt_basal_file, y_snp_file, y_locus_file, y_basal_file, anno_file, bin_choice = boot()

mt_mut_file = "./data/mt_mutation"
mt_basal_file = "./data/mt_phylo"
y_snp_file = "./data/y_snp"
y_locus_file = "./data/y_locus"
y_basal_file = "./data/y_phylo"
anno_file = "./data/anno_file"

# Calling functions
y_df, mt_df = open_data(anno_file)
y_df, mt_df = create_bins(y_df, mt_df, bin_choice)
y_df = findYBasalHaplogroups(y_df, y_snp_file, y_basal_file, y_locus_file)
mt_df = findMTBasalHaplogroups(mt_df, mt_mut_file, mt_basal_file)
y_df, mt_df = createDummyVariables(y_df, mt_df)

# Merge y_df and mt_df into a dataframe carrying the mt- and y-haplogroups for each observation
combined_df_inds = pd.merge(y_df, mt_df, on=['Master ID', 'CombinedBins', 'Lat.', 'Long.'], how='outer')
combined_df_nations, combined_df_inds = createTable(combined_df_inds)
saveToFile(combined_df_inds, combined_df_nations)


# Group by 'Lat.' and 'Long.' and apply custom aggregation
    
################################################
    # Tweaka UI 
    # Ändra readme och artikel
    # Fixa git och github, se till att länka i artikeln
    # 

createWebApplication(combined_df_inds, combined_df_nations)

                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
