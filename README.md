This repository contains the R scripts and necessary data files to estimate the models and run the predictions presented in the paper: "Systematic estimates of the global causes of neonatal and under-five mortality in 2000 - 2024", published in the BMJ in 2026.

How to replicate the predictions:

Follow these steps:

1. PREPARATION of DATA
    1. Download the R-scrips from the main folder of this repository wherever is going to be the working directory and in your local system.
    2. Create the following subfolders in your working directory: DATA, MODELS and PREDICTIONS (keep the big letters)
    3. Download the following data files from the DATA subfolder in this repository to the DATA subfolder in your working directory: "CovariateDatabase2023-wide_20250728_reduced.csv", "CountryModelClass_20240814.xlsx", "CovariateDatabaseIndia2023-wide_20250724.csv", "good_VRneonates_observed_20250828.dta", "data_neonates_20250912.RData".
    4. You can choose to also download into you DATA subfolder the file "data_neonates_LM_20250830.RData", or you can re-create this file in the next step.
    5. (OPTIONAL) If you want to recereate the "data_neonates_LM_20250830.RData" file, just run the script "Prepare data for LM models.R" in your working directory. This script will read the files "good_VRneonates_observed_20250828.dta" and "CovariateDatabase2023-wide_20250728.csv" from the DATA subfolder and will write a file "data_neonates_LM_20250830.RData" into that subfolder. If you want to skip this step, just copy the existing "data_neonates_LM_20250830.RData" from the DATA subfolder in GitHub to your local DATA subfolder.
2. ESTIMATION of MODELS
    1. Generate the Models for LOW Mortality countries: Run the cript "Estimate_LM.R" in your working directory. This will read the file "data_neonates_LM_20250830.RData" from your DATA subfolder and will generate two model files in your MODELS subfolder: "Stan_LM_early_new_0.07_1.RData", "Stan_LM_late_new_0.07_1.RData". Unfortunately these files ar to large (around 62MB each) and we could not provide copies of them here.
    2. Generate the Models for HIGH Mortality countries: Run the cript "Estimate_HM.R" in your working directory. This will read the file "data_neonates_20250912.RData" from your DATA subfolder and will generate two model files in your MODELS subfolder: "Stan_HM_CD_new_unc_0.07_150.RData", "Stan_HM_CD_new_cal_0.07_150.RData". Unfortunately these files ar to large (around 70MB each) and we could not provide copies of them here.
3. PREDICTION of DEATHS
    1. Run the R-script  
