rm(list=ls())
#####################################################################################################################
###
###   Analysis of final data includingseparate Indian states
###
###   J Perin (July 2025)
###
###   Top ten burden countries (1-59 months)
###
### Nigeria (NGA)
### India (IND) ---- envelopes, aids and measles?
### DRC (COD)
### Pakistan (PAK)
### Ethiopia (ETH)
### China (CHN)
### Tanzania (TZA)
### Indonesian (IDN)
### Angola (AGO)
### Bangladesh (BGD)
### Niger (NER)
###
###
#####################################################################################################################

rm(list=ls())
Start<- Sys.time()
options(width=110)




  
#### Analysis of samples
library(reshape2)
library(R2jags)
library(doParallel)
library(coda)
library(lattice)
library(foreign)
library(readstata13)
library(haven)
library(plyr)


# Read data
load("../data/Make_input_data_20250721_calib_JP.RData")

table(studies$nationalrep)
unique(studies$iso3[which(studies$nationalrep==1)])
names(studies)
#studies$iso3[which(studies$nationalrep==1), c("iso3","Refid")]


ls()
studies$rG<- as.numeric(as.factor(studies$iso3))
studies$pcv<- studies$pcv10+studies$pcv13+studies$pcv7
vdr

xt<- which(!is.na(studies$india_state))
summary(studies$meningitis[xt]/studies$totdeaths[xt])


# Whether to use WHO malaria (national history model)
MW<- read.csv("../data/WHOmal_20241009.csv")
table(MW$malwho, useNA="always")
isos.malwho<- unique(MW$iso3[which(MW$malwho==1)])

# remove Namibia because WHO estimate is >60%
isos.malwho<- isos.malwho[which(isos.malwho !="NAM")]

#for malaria and envelopes
path<- "../../other_inputs/"
#list.files(path)
library(foreign)
sca<- read.dta("../../other_inputs/other_inputs_20260224.dta")
summary(sca$igme_pnd-sca$pnd)

sc<- sca[,c("iso3","year","measin","measout","tb_nonresp","tb_resp",
                  "hivp", "mal", "malwho","pnd","u5d")]





# Covarites
path<-"../data/"
nc<- read.dta(paste0(path,"vavr_covariates_20250930.dta"))
summary(nc$sanitation)
table(nc$year)
#nc$pcv
summary(nc$pcv)
table(nc$pcv>1)
# for India states
#nc$iso3[which(nchar(nc$iso3)==2)]
table(is.na(nc$iso3))
nc$iso3<- ifelse(nc$iso3=="NA", nc$india_state23, nc$iso3)
table(nchar(nc$iso3))

# as of 07/21/2025, sc is missing India states,
# so those are dropped from data here
tmp<- merge(sc, nc, by=c("iso3","year"))
unique(tmp$iso3)
summary(tmp$u5mr)
summary(tmp$sanitation)
sc<- tmp

# sc$pnd <- sc$igme_pnd - sc$hivp
sc$envelope_aids_measles_free<- sc$pnd - sc$measin
summary(sc$envelope_aids_measles_free)

head(sc[,c("iso3","year","measin","measout","pnd","u5d","hivp",
                                                 "envelope_aids_measles_free")])

head(sc[which(sc$envelope_aids_measles_free<0),c("iso3","year","measin","measout","pnd","hivp",
                                                 "envelope_aids_measles_free")])
head(sc[which(sc$iso3=="NIU"),c("iso3","year","measin","measout",
                                "pnd","hivp")], n=20)
summary(sc)
# FOR NIUE
sc$envelope_aids_measles_free<- ifelse(sc$envelope_aids_measles_free==0,1,
                                                 sc$envelope_aids_measles_free)

names(sc)
data.predict<- sc
data.predict$isocode<- data.predict$iso3
dim(data.predict)
table(data.predict$year)
length(unique(data.predict$isocode))
summary(data.predict$envelope_aids_measles_free==0)


###################################X
#######   PREDICTIONS     ##########
###################################X
load(file=paste0("../results/VA_lambda400_rsdlim0.07_20250731.RData"))
ebn<-colMeans(out$MCMCout$B)
ebn

summary(out$studies$pcv)
summary(data.predict$pcv)

summary(out$studies$rota)
summary(data.predict$rota)

vxf<- out$vxf
vxf
vdr

studies$intercept<-1
data.predict$intercept<-1

for(i in vxf){
  print(i)
  print(summary(studies[,i]))
  print(summary(data.predict[,i]))
}




#function(MO, STUP, NP=500, IDV=c("isocode","year"), PPR=c(0, 0))
#iterations saved ,  chains, parameters
dim(out$MCMCout$B)

# Scale covariates to match scale used in studies
XX <- data.predict[, c(vxf)] 
for(vx in vxf[3:length(vxf)]){
  mn<- out$Vars$mean[which(out$Vars$xvar==vx)]
    stddev<-out$Vars$sd[which(out$Vars$xvar==vx)]
  XX[,vx]<- (XX[, vx] - mn)/stddev
}


# Random effects
# Assign random effect to country years
# only use if there is nationally representative study 
MR <- colMeans(out$MCMCout$re, na.rm = T)    # unweighted means of RE
rownames(MR)
iso.natrep<- unique(studies$iso3[which(studies$nationalrep==1)])
iso.natrep<- iso.natrep[which(iso.natrep!="NER")]
names(studies)
studies$iso3

findRE<- function(iso) { 
#  print(iso)
  if(iso %in% iso.natrep) re<- MR[iso,]
  else re<- rep(0, length(vdr)-1)
  re 
}
findRE("AFG")
RE<- do.call("rbind", lapply(data.predict$iso3, findRE))


# Determine fixed effect and include random effects
# Multiply every iteration as matrix with prediction data
FFa <- array(NA, dim=c(nrow(data.predict), length(vdr), nrow(out$MCMCout$B)) )
FFa[] <- apply(out$MCMCout$B, 1, function(x) {
  M<-exp(cbind(0, as.matrix(XX) %*% as.matrix(x) + RE) )
  # translate linear predictors to probabilities
  MP<- M / rowSums(M)
  return(MP) }
  )
dim(FFa)

rownames(FFa)<- paste0(data.predict$iso3,".", data.predict$year)
colnames(FFa)<- vdr
head(FFa[,,1])
# table(rowSums(FFa[,,1]))




####################################
#######   WHO REGIONS   ############
####################################
codes<- read.dta13("../../other_inputs/codelist region.dta")
#head(codes)
codes$isocode<- codes$iso3


####################################
#######   WHO MALARIA   ############
####################################
table(data.predict$iso3)
data.predict$iso.nat<- ifelse(nchar(data.predict$isocode)>3, "IND", data.predict$isocode)
unique(data.predict$iso3)
unique(data.predict$iso.nat)



library(arrayhelpers)
NSim<- dim(FFa)[3]

############################################################
rownames(FFa)[1:20]
xt<- paste0("NGA.",c(2000:2023))
# xt<- paste0("PSE.",c(2000:2023))
# FFa[xt,,1]

xt<- which(data.predict$pfpr==0)
summary(FFa[xt,,1])





############################################################
testIND<- FFa

summary(data.predict$envelope_aids_measles_free)
dim(testIND)
dimnames(testIND)

xt<- paste0("AGO.",c(2000:2023))
testIND[xt,,1]



############################################################
testA<- testIND

#Make data set of country point estimates
dim(testIND)
est<-apply(testA, c(1,2), FUN=mean)
#lower<-apply(testA, c(1,2), FUN=quantile, probs=c(0.025))
#upper<-apply(testA, c(1,2), FUN=quantile, probs=c(0.975))
colnames(est) <- vdr #colnames(lower)<-colnames(upper)<-vdr
dim(est)
dim(data.predict)
table(data.predict$iso3)
rownames(est) <- paste0(data.predict$isocode,".", data.predict$year) #rownames(p21$estR)
dim(est)
head(est)
VA.est<- as.data.frame(est) #merge(est, lower, by="row.names")
head(VA.est, n=30)

temp<-strsplit(rownames(VA.est), split="\\.")
head(temp)
iso<-as.vector(unlist(lapply(temp , FUN= function(x) x[[1]])))
year<-as.vector(unlist(lapply(temp , FUN= function(x) x[[2]])))
VA.est$iso3<- iso #dp.nat$isocode
VA.est$year <- year #dp.nat$year
dim(VA.est)
VA.est[which(VA.est$iso3=="IDN"),]
VA.est[which(VA.est$iso3=="IND"),]
VA.est[which(VA.est$iso3=="IOri"),] # check India states



dim(VA.est)
table(VA.est$iso3)


#VA.est<- merge(VA.est, codes[,c("iso3","whoreg6")], by="iso3")
head(VA.est)
dim(VA.est)

VA.est$Row.names<- NULL
head(VA.est)
VA.draws<- testA
VA.draws.noPH<- FFa

#test<-cbind(VA.est$iso3, dp.nat$isocode)
#test[900:1200,]

#fractions
VA.est[which(VA.est$iso3=="NGA"),]

head(VA.est)






table(VA.est$year)
dim(VA.est)



#COEFFICIENTS
np<- length(vxf) 
Best<- colMeans(out$MCMCout$B)
Best

length(rep(2:length(vdt), each=np))
coef<- data.frame( cause = rep(vdt[-1],each=np), cause.n=rep(2:length(vdt), each=np),
                   cov=rep(c(vxf), length(vdt) -1),
                   est=as.vector(Best))
coef

############################################
dim(VA.est)
save(VA.draws, VA.est,
     file=paste0("../results/Predict_VA.RData"))
