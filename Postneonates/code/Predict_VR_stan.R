
########################################
###
###   Analysis of final data including two periods for separate indian states
###
###   David Prieto (Feb 2020)
###
###
###
########################################


rm(list=ls())
Start<- Sys.time()
options(width=110)


#### Analysis of samples
library(R2jags)
library(doParallel)
library(coda)
library(lattice)
library(foreign)
library(readstata13)
#library(saveJAGS)
library(readstata13)

###
###  Prepare datasets for R
###

# Read data
# new March 2025
load("../data/Make_input_data_VR.RData") 

#vxf #covariates
vdt #causes
dim(studies)
table(studies$year)
data.predict<- read.dta13("../data/vavr_covariates_20250930.dta")
table(data.predict$year)

#for envelopes, TB, measles, malaria
path<- "../../other_inputs/"
#list.files(path)
library(foreign)
sca<- read.dta(paste0(path,"other_inputs_20260224.dta"))

sc<- sca[,c("iso3","year","measin","measout","tb_nonresp","tb_resp",
            "hivp", "mal", "malwho","pnd","u5d","igme_pnd")]
sc$envelope_aids_measles_free<- sc$pnd - sc$measin

# Need to add India to other_inputs
xt<- merge(data.predict, sc, by=c("iso3","year"))
table(data.predict$year)
table(xt$year)
data.predict<- xt
summary(data.predict$envelope_aids_measles_free)

# for NIUE 2005
data.predict$envelope_aids_measles_free<- ifelse(data.predict$envelope_aids_measles_free==0,1,
                                                 data.predict$envelope_aids_measles_free)
dim(data.predict)
table(data.predict$year)



###################################X
#######   PREDICTIONS     ##########
###################################X

#Contains object named "out"
#load(file="../results/VR_Bayes_lambda5_sd_0.07.RData")
# load(file="../results/VR.RData")
#load(file="../results/VR_lambda0.9_rsdlim0.04.RData")
load(file="../results/VR_lambda0.9_rsdlim0.04.RData")
ebn<-colMeans(out$MCMCout$B)
ebn

vxf<- out$vxf
vxf
vdr

studies$intercept<-1
data.predict$intercept<-1
data.predict$Iwhoreg6_3<- data.predict$Iwhoreg6_SEAR
#data.predict$Hib3_cov<- data.predict$vac_hib3

for(i in vxf){
  print(i)
  print(summary(studies[,i]))
  print(summary(data.predict[,i]))
}

table(data.predict$year)


# Random effects
# Assign random effect to country years
# only use if there is nationally representative study 
MR <- colMeans(out$MCMCout$re, na.rm = T)    # unweighted means of RE
rownames(MR)
head(MR)
summary(MR)
dim(out$MCMCout$re)
head(out$MCMCout$re[1,,])
head(out$MCMCout$re[1,2,])



#function(MO, STUP, NP=500, IDV=c("isocode","year"), PPR=c(0, 0))
#iterations saved ,  chains, parameters
dim(out$MCMCout$re)
head(out$MCMCout$re)
dim(MR)
nrow(out$MCMCout$re)
# assign random RE, b/c
# countries with data are not modeled
# need to have RE to get best variance
findRE<- function(iso) { 
  #  print(iso)
  # select random iteration
    it<- sample(1:nrow(out$MCMCout$re),1)
  # select random country
    rn<- sample(1:nrow(MR),1)
    re<- out$MCMCout$re[it, rn,]
  re 
}
RE<- do.call("rbind", lapply(data.predict$iso3, findRE))


# Scale covariates to match scale used in studies
XX <- data.predict[, c(vxf)] 
for(vx in vxf[3:length(vxf)]){
  mn<- out$Vars$mean[which(out$Vars$xvar==vx)]
  stddev<-out$Vars$sd[which(out$Vars$xvar==vx)]
  XX[,vx]<- (XX[, vx] - mn)/stddev
}


# Determine fixed effect and include random effects
# Multiply every iteration as matrix with prediction data

FFa <- array(NA, dim=c(nrow(data.predict), length(vdr), nrow(out$MCMCout$B)) )
FFa[] <- apply(out$MCMCout$B, 1, function(x) {
  REs<- do.call("rbind", lapply(data.predict$iso3, findRE))
  M<-exp(cbind(0, as.matrix(XX) %*% as.matrix(x) + REs) )
  # translate linear predictors to probabilities
  MP<- M / rowSums(M)
  return(MP) }
)
dim(FFa)

rownames(FFa)<- paste0(data.predict$iso3,".", data.predict$year)
colnames(FFa)<- vdr
head(FFa[,,1])
dim(FFa)
# table(rowSums(FFa[,,1]))







library(arrayhelpers)
library(plyr)

dim(FFa)
rownames(FFa)
head(FFa[grepl("SAU",rownames(FFa)),,1], n=25)
est<-apply(FFa, c(1,2), FUN=mean, na.rm=TRUE)
dim(est)
#lower<-apply(p21$estR, c(1,2), FUN=quantile, probs=c(0.025), na.rm=TRUE)
#upper<-apply(p21$estR, c(1,2), FUN=quantile, probs=c(0.975), na.rm=TRUE)
colnames(est) <- vdr #colnames(lower)<-colnames(upper)<-vdr
#rownames(est) <-rownames(p21$estR)
#colnames(lower)<- paste(colnames(lower),".lower",sep="")
#colnames(upper)<- paste(colnames(upper),".upper",sep="")

#dim(est)
#colnames(est)
test<- rowSums(est)
table(test)

VR.est<- as.data.frame(est)
#rownames(VR.est)<- VR.est$Row.names
temp<-strsplit(rownames(VR.est), split="\\.")
iso<- unlist(lapply(temp, FUN=function(x) x[[1]]))
year<- unlist(lapply(temp, FUN=function(x) x[[2]]))
VR.est$iso3<- iso #data.predict$isocode
VR.est$year <- year #data.predict$year
VR.est$whoreg6 <- data.predict$whoreg6
head(VR.est)

VR.est[which(VR.est$iso3=="IDN"),]

head(VR.est)
VR.draws<- FFa #p21$estR


summary(VR.est)
table(VR.est$year)
dim(VR.est)
dim(VR.draws)
rownames(VR.draws)
save(VR.draws, VR.est, file=paste0("../results/Predict_VR.RData"))
table(VR.est$year)

head(data.predict[which(VR.est$iso3=="SAU"),vxf])

#COEFFICIENTS
#out$MCMCout$B
data.predict[which(data.predict$iso3=="IDN"), c("iso3","year",vxf)]
colMeans(out$MCMCout$B)

head(VR.est[which(VR.est$iso3=="SAU"),])
np<- length(vxf)
coef<- data.frame( cause = rep(vdt[-1],each=np), cause.n=rep(2:8, each=np),
                   cov=rep(c(vxf), 7), est=as.vector(colMeans(out$MCMCout$B)))
head(coef)
if(!interactive()) q("no")
vxf
vdt
