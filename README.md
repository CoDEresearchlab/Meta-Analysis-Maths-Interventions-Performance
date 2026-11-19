# Meta-Analysis-Maths-Interventions-Performance

This meta-analysis examines the extent to which interventions improve mathematics performance in primary education. We synthesised evidence from 280 studies that used randomised controlled trials or quasi-experimental designs, analysing a total of 614 effect sizes from 82,047 participants. 

### Repository contents

This repository is organised into one folder only. The folder contains the data and all R scripts used to calculate effect sizes and generate statistics, figures, and tables.

### Data source

#### Raw data
The raw data are stored in `mydata.xlsx`. This file contains all the information needed to calculate effect sizes for each study.

#### Processed data
After calculating effect sizes, the processed data are stored in `mydata_processed.csv`.

#### Data for running meta-regression models
Rows with insufficient data for running the models are removed, and the resulting dataset is stored in `mydata_cog.csv`.

### R scripts

In general, R scripts are named according to the figures, tables, or analyses they produce.

#### 1. Effect size calculation.R
This script calculates the effect sizes for each study using the raw data and generates the processed dataset.

#### 2. Figure 2, S1, S2.R
This script produces Figure 2, Supplementary Figures S1 and S2, and the associated statistical results.

#### 3. Figure 3, suppl table S4, S5.R
This script produces Figure 3, Supplementary Tables S4 and S5, and the associated statistical results.

#### 4. Figure 4.R
This script produces Figure 4 and the associated statistical results.

#### 5. Figure 5, S9.R
This script produces Figure 5, Supplementary Figure S9, and the associated statistical results.

#### 6. Suppl table S7 & results of effectiveness of features.R
This script produces Supplementary Table S7 and the statistical results related to the effectiveness of intervention features.

#### 7. Suppl figure S3, S5-S7 & suppl notes 2-5.R
This script produces Supplementary Figures S3 and S5–S7 and the analyses reported in Supplementary Notes 2–5.

#### 8. Suppl figure S4.R
This script produces Supplementary Figure S4 and the associated statistical results.

#### 9. Suppl table 11.R
This script produces Supplementary Table S11 and the associated statistical results.

#### 10. Suppl table S1, S3.R
This script produces Supplementary Tables S1 and S3 and the associated statistical results.

#### 11. Suppl table S6.R
This script produces Supplementary Table S6 and the associated statistical results.

#### 12. Suppl table S10.R
This script produces Supplementary Table S10 and the associated statistical results.
