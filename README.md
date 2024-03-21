# HaploMapper v1.0

HaploMapper provides a graphical interface for you to explore the distribution of basal haplogroups present in the AADR dataset ([https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data)).

## HaploMapper Structure

- `boot()`: Loads the files present in the data folder.
- `create_bins()`: Bins the AADR dataset based on political entity (country) and sample age. After binning with respect to both columns, a combinedBin is generated (referred to as CombinedBins).
- `FindYHaplogroups()` and `FindMTHaplogroups()`: Finds the basal haplogroups of the ones provided in the AADR dataset. This is achieved by consulting the y_snp, y_locus, y_phylo, and mt_mut, mt_phylo for yDNA and mtDNA respectively (all files are present in the data folder).
- `createDummyVariables()`: Ensures all haplogroups are alphanumerical, then creates a binary representation of the basal haplogroups and sums the occurrences for each sample.
- `createTable()`: Bins the samples with respect to CombinedBins and generates frequency tables for yDNA and mtDNA. Each row in the dataframe outputted from this function corresponds to the basal haplogroups distribution in a political entity (country) during a specific time interval.
- `saveToFile()`: Generates files corresponding to the frequency tables in CSV format.
- `createWebApplication()`: Uses Dash to generate an HTML dashboard that launches the graphical interface of HaploMapper.

## Usage

HaploMapper is run as follows: `python HaploMapper.py`

HaploMapper is coded in Python 3.10.12

## Dependencies

To install HaploMapper's dependencies, you can create a `requirements.txt` file containing the following lines, and then install all dependencies at once using `pip install -r requirements.txt`:

