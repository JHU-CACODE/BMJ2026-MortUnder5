# Estimate high mortality models for neonates
# David Prieto-Merino
# September 2025


rm(list=ls())
library(readxl)
library(tidyverse)
library(flextable)
library(rmdformats)
library(confintr)
library(MCMCvis)
library(bayesplot)
library(rstudioapi)
library(callr)
library(Hmisc)
library(ggpubr)
library(sp)
library(rstatix)
require(rstan)




# Load data prepared previously with the code "Prepare data or models.R" 
load("DATA/data_neonates_20250912.RData")



# Define function for estimations

estimate_HM <- function(STUD, DEAT, SDRE, LAMB, NAME, PATH="MODELS"){
  # STUD: study data.frame
  # DEAT: deaths data.frame
  # SDRE: Max SD of random effects (hyperparameter)
  # LAMB: Lambda coef of Lasso models (hyperparameter)  
  # NAME: Prefix for name of model e.g. "Stan_LM_late_"
  # create input list in global environment
  .GlobalEnv$st.input <- list(
    model = stan_model(file='cacode_DP.stan', auto_write = TRUE),
    studies = arrange(STUD, iso3, recnr) %>%
      mutate(iso3=factor(iso3, levels=sort(unique(iso3))), intercept=1) %>%
      dplyr::select(all_of(c("recnr","strata_id", "iso3", "first", "last", "intercept", vxc, vxn))) %>%
      as_tibble(),
    deaths = dplyr::select(DEAT, all_of(c("row","recnr","cause","n",vdt))) %>% as_tibble(),
    lambda = LAMB,  # lambda for lasso
    rsdlim = SDRE,  # define max SD of RE 
    sdbeta = 0.5,   # SD of parameters not in Lasso
    vdt = vdt,
    vxn = vxn,
    vxc = c("intercept",vxc),
    nsv = 3,
    param = c('B', 're', 'sd_re'),
    nchai = 4,
    niter = 4000,
    nwarm = 2000,
    cores = 4,
    patho = PATH, # path to store model output
    name = paste(NAME,SDRE,LAMB, sep="_"), # name of the model
    summary = 0  # put 0 if you want to save all the posterior samples of the parameters to estimate, or 1 if you only want to save the summaries of the posterior distributions  
  )
  # Run job on background 
  jobRunScript("runStan2.R", importEnv=TRUE, exportEnv="R_GlobalEnv")
}



estimate_HM(studies_23, mutate(deaths_23, n=unc), 0.07, 150, "Stan_HM_CD_new_unc")
estimate_HM(studies_23, mutate(deaths_23, n=cal), 0.07, 150, "Stan_HM_CD_new_cal")

# save(Stan_HM_CD_new_cal_0.07_150, file="../results/Stan_HM_CD_new_cal_0.07_150.RData")








