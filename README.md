<!--your zenodo badge here-->

# wild-etal_2022_tbd
Meta-Repository for Publication: \
**Large EMulated ENSembles (LEMENS) Enable Robust Evauation of Climate Impacts on Extreme Droughts Globally**

[![DOI](https://zenodo.org/badge/265119113.svg)](https://zenodo.org/badge/latestdoi/265119113)

Thomas B. Wild<sup>1,2,3\*</sup>, Claudia Tebalid<sup>1</sup>, Chris Vernon<sup>1</sup>, Abigail Snyder<sup>1</sup>, Kalyn Dorheim<sup>1</sup>, Mohamad Hejazi<sup>4</sup>, Kamal Chowdhury<sup>2</sup>

<sup>1 </sup> Joint Global Change Research Institute, Pacific Northwest National Laboratory, College Park, MD, USA
<sup>2 </sup> Earth System Science Interdisciplinary Center (ESSIC), University of Maryland, College Park, MD, USA
<sup>3 </sup> Department of Civil and Environmental Engineering, University of Maryland, College Park, MD USA
<sup>4 </sup> King Abdullah Petroleum Studies and Research Center, Riyadh, Saudi Arabia

\* corresponding author:  thomas.wild@pnnl.gov

## Purpose
This meta-repository creates a single point of access for interested researchers to access all of the components that were used to the published work referenced above for the purpose of reproducibility. This repository contains references to all minted data and software and houses ancillary code used to transform the source data, create figures for the publication, conduct the experiment, and execute the contributing software.

## Abstract
In this paper the authors develop a long-term global energy-economic model which is capable of assessing alternative energy evolutions over periods of up to 100 years. The authors have sought to construct the model so that it can perform its assigned task with as simple a modelling system as possible. The model structure is fully documented and a brief summary of results is given.

## Journal reference
Edmonds, J., & Reilly, J. (1983). A long-term global energy-economic model of carbon dioxide release from fossil fuel use. Energy Economics, 5(2), 74-88. DOI: https://doi.org/10.1016/0140-9883(83)90014-2

## Code reference
References for each minted software release for all code involved.  

These are generated by Zenodo automatically when conducting a release when Zenodo has been linked to your GitHub repository. The Zenodo references are built by setting the author order in order of contribution to the code using the author's GitHub user name.  This citation can, and likely should, be edited without altering the DOI.

If you have modified a codebase that is outside of a formal release, and the modifications are not planned on being merged back into a version, fork the parent repository and add a `.<shortname>` to the version number of the parent and construct your own name.  For example, `v1.2.5.hydro`.

Human, I.M. (2021, April 14). Project/repo:v0.1.0 (Version v0.1.0). Zenodo. http://doi.org/some-doi-number/zenodo.7777777

## Data reference

### Input data
Reference for each minted data source for your input data.  For example:

Human, I.M. (2021). My input dataset name [Data set]. DataHub. https://doi.org/some-doi-number

### Output data
Reference for each minted data source for your output data.  For example:

Human, I.M. (2021). My output dataset name [Data set]. DataHub. https://doi.org/some-doi-number

## Contributing modeling software
| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| model 1 | version | link to code repository | link to DOI dataset |
| model 2 | version | link to code repository | link to DOI dataset |
| component 1 | version | link to code repository | link to DOI dataset |

## Reproduce my experiment
Fill in detailed info here or link to other documentation that is a thorough walkthrough of how to use what is in this repository to reproduce your experiment.

1. Install the software components required to conduct the experiement from [Contributing modeling software](#contributing-modeling-software)
2. Download and install the supporting input data required to conduct the experiement from [Input data](#input-data)
3. Run the following scripts in the `workflow` directory to re-create this experiment:

| Script Name | Description | How to Run |
| --- | --- | --- |
| `step_one.py` | Script to run the first part of my experiment | `python3 step_one.py -f /path/to/inputdata/file_one.csv` |
| `step_two.py` | Script to run the last part of my experiment | `python3 step_two.py -o /path/to/my/outputdir` |

4. Download and unzip the output data from my experiment [Output data](#output-data)
5. Run the following scripts in the `workflow` directory to compare my outputs to those from the publication

| Script Name | Description | How to Run |
| --- | --- | --- |
| `compare.py` | Script to compare my outputs to the original | `python3 compare.py --orig /path/to/original/data.csv --new /path/to/new/data.csv` |

## Reproduce my figures
Use the scripts found in the `figures` directory to reproduce the figures used in this publication.

| Script Name | Description | How to Run |
| --- | --- | --- |
| `generate_figures.py` | Script to generate my figures | `python3 generate_figures.py -i /path/to/inputs -o /path/to/outuptdir` |
