# fNIRS validation study

This repository contains files for producing the results found in: (to be added)

The data repository for this study can be found in: https://doi.org/10.48723/vscr-eq07. Details on how this dataset was prepared can be found in: https://github.com/alkvi/fnirs_dataset_preparation.

## Data analysis

Analysis is performed on the data found in the data repository.

These scripts assume the following directory structure:

fnirs_validation_study  
│--   Protocol_1_analysis.m  
│--   ...  
Data  
│--   (fNIRS data etc)  
│--   temp_data  
│--   mat_files  

The temp_data and mat_files folders have to be created. Note that the BIDS dataset comes zipped and needs to be unzipped before analysis.

### Demographics and gait variables

Important demographics are created and then gait variables used in the mixed models are extracted from the full dataset.

~~~
python summarize_demographics.py
python prepare_gait_params_for_mixed_model.py
~~~

### Signal quality

To run the quality scripts, run them from the signal quality directory
~~~
cd Signal_quality
~~~

This calculates SCI per condition and participant, and calculates channels below acceptable threshold (0.7), and bad channels are formatted for analysis.
~~~
python signal_quality_complex_walk.py
python prepare_bad_channel_list.py
~~~

Peak power is calculated by running _power_calculation.m_ in MATLAB. 

### 1st and 2nd level analysis

These steps are all run in MATLAB.

- Preprocessing and 1st level analysis is performed in Protocol_1_analysis.m
- 2nd level analysis is performed in Protocol_1_group_analysis.m. This is run on each hemoglobin type, specified in the beginning of the script.

## Producing figures

- A figure with visualized T statistics for each hypothesis is created with the R markdown file _summarize_results.Rmd_.
- To prepare 1st level ROI beta variables for visualization, the file _prepare_betas_interaction_plot.m_ is run.
- A figure with 1st level beta variables plotted against gait variables is created with _plot_interactions.Rmd_.

