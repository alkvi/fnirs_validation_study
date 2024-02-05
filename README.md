# fNIRS validation study

This repository contains files for producing the results found in: (to be added)

The data repository for this study can be found in: https://doi.org/10.48723/vscr-eq07. Details on how this dataset was prepared can be found in: (to be added).

## Data analysis

Analysis is performed on the data found in the data repository.

### Demographics and gait variables

Important demographics are created using _summarize_demographics.py_.

Gait variables used in the mixed models are extracted from the full dataset in _prepare_gait_params_for_mixed_model.py_.

### Signal quality

- Signal quality in terms of SCI is evaluated using _Signal_quality/signal_quality_complex_walk.py_. This calculates SCI per condition and participant, and calculates channels below acceptable threshold (0.7).
- Peak power is calculated using _Signal_quality/power_calculation.m_.
- Bad channels are formatted for analysis using _Signal_quality/prepare_bad_channel_list.py_.

### 1st and 2nd level analysis

- Preprocessing and 1st level analysis is performed in Protocol_1_analysis.m
- 2nd level analysis is performed in Protocol_1_group_analysis.m. This is run on each hemoglobin type, specified in the beginning of the script.

## Producing figures

- A figure with visualized T statistics for each hypothesis is created with the R markdown file _summarize_results.Rmd_.
- To prepare 1st level ROI beta variables for visualization, the file _prepare_betas_interaction_plot.m_ is run.
- A figure with 1st level beta variables plotted against gait variables is created with _plot_interactions.Rmd_.

