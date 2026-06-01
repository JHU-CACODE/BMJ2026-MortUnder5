# ---
#   title: "Estimation with Stan models of LM models"
# author: "David Prieto-Merino"
# date: "`r Sys.time()`"
# output:
#   html_document:
#   toc: true
# df_print: paged
# self.contained: true
# word_document:
#   toc: true
# pdf_document:
#   toc: true
# highlight: tango
# includes:
#   toc_depth: 4
# editor_options:
#   chunk_output_type: console
# ---
  

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


# Load NEW data prepared previously with the code "Prepare data or models.R" 
load("DATA/data_neonates_LM_20250830.RData")


# Define function for estimations

estimate_LM <- function(STUD, DEAT, SDRE, LAMB, NAME, PATH="MODELS"){
  # STUD: study data.frame
  # DEAT: deaths data.frame
  # SDRE: Max SD of random effects (hyperparameter)
  # LAMB: Lambda coef of Lasso models (hyperparameter)  
  # NAME: Prefix for name of model e.g. "Stan_LM_late_"
  
  # create input list in global environment
  .GlobalEnv$st.input <- list(
    model = stan_model(file='cacode_DP.stan', auto_write = TRUE),
    studies = arrange(STUD, iso3, recnr) %>%
      mutate(iso3=factor(iso3, levels=sort(unique(iso3)))) %>%
      dplyr::select(all_of(c("recnr", "country", "iso3", "first", "last", "intercept", vxn))) %>%
      as_tibble(),
    deaths = dplyr::select(DEAT, all_of(c("row","recnr","cause","n",vdt))) %>% as_tibble(),
    lambda = LAMB,  # lambda for lasso
    rsdlim = SDRE,  # define max SD of RE 
    sdbeta = 0.5,   # SD of parameters not in Lasso
    vdt = vdt,
    vxn = vxn,
    vxc = "intercept",
    nsv = 1,
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


estimate_LM(studies_e,  deaths_e,  0.07, 1, "Stan_LM_early_new")
estimate_LM(studies_l,  deaths_l,  0.07, 1, "Stan_LM_late_new")


# save(Stan_LM_early_new_0.07_1, file="../results/Stan_LM_early_new_0.07_1.RData")
# save(Stan_LM_late_new_0.07_1, file="../results/Stan_LM_late_new_0.07_1.RData")







