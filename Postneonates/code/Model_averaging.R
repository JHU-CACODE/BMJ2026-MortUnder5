
########################################
###
###   Analysis of final data including two periods for separate indian states
###
###   David Prieto (Feb 2020)
###
########################################

rm(list=ls())
library(coda)
library(lattice)
library(foreign)
library(tidyverse)
library(ggpubr)
library(abind)
library(readstata13)


type_c<- ""


#Predictions from VA model
# wasting, hib, rota
#load("../results/Predict_VA_20250717.RData")
load(paste0("../results/Predict_VA.RData"))
length(unique(VA.est$iso3))
head(VA.est)
table(VA.est$year)

#Predictions from VR model
load("../results/Predict_VR.RData")
head(VR.est[which(VR.est$iso3=="SAU"),])
length(unique(VR.est$iso3))

names(VR.est)<- gsub("perinatal","neonatal", names(VR.est))
head(VR.est)
table(VR.est$year)

unique(VA.est$iso3) [!unique(VA.est$iso3) %in% unique(VR.est$iso3)]


# load a VA and VR models to get variables needed for prediction
# VR model
#load("../results/VR.RData")
load(file="../results/VR_lambda0.9_rsdlim0.04.RData")
vr.out<- out
vr.out$vxf
vr.out$vdt
it<- rownames(VA.draws)
table(nchar(rownames(VA.draws)))
xt<-nchar(rownames(VA.draws))
it<- substr(it, 1,ifelse(xt==8,3,4))

# VA model
#load("../results/VA_20250526_50.RData")
load(paste0("../results/VA_lambda400_rsdlim0.07_20250731.RData"))
va.out<- out
va.out$vxf
va.out$vdt


## Load data from all countries:
# dc <- read.dta13(file="../data/envelopes_IGME_20220928.dta")
# dc <- read.dta13(file="../data/master national level data_20220928_1980-2020.dta")
# dc <- read.dta13(file="../data/vavr_national_covariates_2022 11 06.dta")


#dc<- read.dta13(file="../data/model_averaging_weight_2023 02 14.dta")
dc<- read.dta13(file="../../other_inputs/other_inputs_20260224.dta")
length(unique(dc$iso3))
names(dc)
dc<- dc[,c("iso3","year","MA_weight")]
head(dc[which(dc$iso3=="SAU"),])
#dc$u5mr<- dc[,"_5q0"]
summary(dc$igme_u5mr)
table(dc$iso3[which(is.na(dc$igme_u5mr))])
table(dc$iso3)
table(dc$year)
length(unique(it))
uit<- unique(it)
udc<- unique(dc$iso3)
table(uit %in% udc)
uit[ !(uit %in% udc) ]

dp<- dc[which(dc$iso3 %in% it & dc$year>=2000),]
names(dp)
## di <- read_csv(file="india_covs_1990-2019_1May2020.csv")

###  Combine VA and VR predictions
###
## For each country/year crate a weight of VA dep
table(dp$iso3)
dp$isoyear <- paste(dp$iso3, dp$year, sep=".")
dp$wei.vavr <- sapply(dp$MA_weight, function(x) max(0, min(1, (x - 25)/ (35 - 25) )))
dp[which(dp$iso3=="IDN"), c("iso3","year","MA_weight")] #, "ln_5q0_cov")]
summary(dp$wei.vavr)

# causes in each model
ca <- va.out$vdt

cr <- vr.out$vdt
cr<- gsub("perinatal","neonatal", cr)
cc <- unique(c(cr,ca))
cc

table(VA.est$year)

## Complete matrices with non-existing causes
# in VA matrices
dim(VA.draws)
rownames(VA.draws[,,1])
cc
ca
PFa <- abind(VA.draws, array(0, dim=c(dim(VA.draws)[1],
                                      length(cc[!(cc %in% ca)]), dim(VA.draws)[3])), along=2)
dim(PFa)
dimnames(PFa)[[2]] <- c(ca, cc[!(cc %in% ca)])
dimnames(PFa)[[1]]
#dp<- dp[which(dp$iso3!="IND"),]
PFa <- PFa[,cc,]
PFa <- PFa[dp$isoyear,cc,]
dim(PFa)

table(rownames(PFa) %in% dp$isoyear)
table( dp$isoyear %in% rownames(PFa) )
dp$isoyear[! dp$isoyear %in% rownames(PFa) ]
table(dp$iso3)
table(dp$year)
rownames(PFa)
PFa[1:3,,1:2]
dim(VA.draws)
dim(PFa)

length(unique(dp$iso3))
PRa <- abind(VA.draws, array(0, dim=c(dim(VA.draws)[1],
                                      length(cc[!(cc %in% ca)]), dim(VA.draws)[3])), along=2)
dimnames(PRa)[[2]] <- c(ca, cc[!(cc %in% ca)])

PRa <- PRa[dp$isoyear,cc,]
PRa[1:3,,1:2]


# in VR matrices
# placeholder - VR model not run yet
#dim(VR.draws)
#VR.draws<- VA.draws
dim(VR.draws)
cc
cr
PFr <- abind(VR.draws, array(0, dim=c(dim(VR.draws)[1],
                                      length(cc[!(cc %in% cr)]),
                                      dim(VR.draws)[3])), along=2)
colnames(PFr)
dim(PFr)
colnames(PRa)
dim(PRa)
#dim(VA.draws)
#dimnames(VA.draws)[[2]]
#dimnames(PFr)[[2]]
dimnames(PFr)[[2]] <- c(cr, cc[!(cc %in% cr)])
PFr <- PFr[dp$isoyear,cc,]
#dimnames(PFr)[[2]]
names(dp)
#table(dp$isoyear %in% dimnames(PFr)[[1]])
cc
#table( dp$isoyear %in% dimnames(PFr)[[1]] )
#dp$isoyear[which ( !dp$isoyear %in% dimnames(PFr)[[1]] )]



#PFr[1:3,,1:2]
PRr <- abind(VR.draws, array(0, dim=c(dim(VR.draws)[1],
                                      length(cc[!(cc %in% cr)]),
                                      dim(VR.draws)[3])), along=2)
dim(PRr)
dimnames(PRr)[[2]] <- c(cr, cc[!(cc %in% cr)])
PRr <- PRr[dp$isoyear,cc,]


## array of weights
WM <- array(dp$wei.vavr, dim=dim(PFa))
head(WM)


xt<-dimnames(PFa[which(!WM[,1,5] %in% c(0,1)),,1])[[1]]

## Combination
dim(PFa)
dim(PFr)
dim(WM)


it<-dimnames(PFr)[[1]]
it2<-dimnames(PFa)[[1]]
it[!it %in% it2]

PFr<- PFr[it2, ,]

pco.a <- list(Model="VAVR",
              PF=(PFa*WM + PFr*(1-WM)), PR=(PRa*WM + PRr*(1-WM)))
MA.draws<- pco.a
save(MA.draws, file="../results/PN.MA.draws.RData")
# Calculate predictions
est<-apply(pco.a$PF, c(1,2), FUN=mean)


dim(est)
colnames(est)
summary(est)

est<- as.data.frame(est)
rownames(est)
colnames(est)
table(nchar(rownames(est)))
xt<- nchar(rownames(est))
est$iso3<- substr(rownames(est),1,ifelse(xt==9,4,3))
est$year<- substr(rownames(est),ifelse(xt==9,6,5),xt)
table(est$year)
#rownames(est)
save(est,  MA.draws, 
     file=paste0("../results/Model_averaging_draws.RData"))
#write.csv(est, file="Prediction2_model_averaging_JP.csv")
#write.dta(as.data.frame(est), file="../results/Prediction2_model_averaging_JP.dta")
table(est$year)
dim(est)

if(!interactive()) q('no')
dim(est)
rownames(est)
summary(est)
head(est, n=25)
dim(est)
rownames(est)

#unique country years
temp<-strsplit(rownames(est), split="\\.")
iso<- unlist(lapply(temp, FUN=function(x) x[[1]]))
year<- unlist(lapply(temp, FUN=function(x) x[[2]]))
unique(iso)
unique(year)
table(year)
table(iso, year==2024)
est[which(grepl("NAM", rownames(est))),]
