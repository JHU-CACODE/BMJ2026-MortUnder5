
options(width=110)
options("width")
#head(studies)
sessionInfo()
rm(list=ls())

## real data ----
### load ----
# rm(list=ls())
load("../data/Make_input_data_VR.RData") 
dim(studies)

table(studies$year)
unique(studies$isocode)

#dim(data.predict)
studies$rG<- as.numeric(as.factor(studies$isocode))
dim(studies)
table(studies$Iwhoreg6_3)
studies$intercept<-1
#load("../results/VR.RData")
#out$Studies$rG
#LASSO_VR_jags = out



# ######################################################################## #
# user specifications ----
# ######################################################################## #

## paths ----
### sourcecode ----
#sourcecode.path = "/Users/sandipanpramanik/mystaff/JHU/research/postdoc/ca-code/cacode modeling scripts/code/sourcecode/5-12-25"

### data ----
#oldcode.path.jags = "/Users/sandipanpramanik/mystaff/JHU/research/postdoc/ca-code/cacode modeling scripts/postneonatal_cod/results"

### where to save output ----
output.path = file.path(getwd(), "../results")

### causes ----
(causes = c("pneumonia", "diarrhea", "meningitis", "perinatal",
            "othergroup1", "congenital", "otherncd", "injuries"))

### covariates ----
(covariates = c("intercept","Iwhoreg6_3", "year","ln_5q0_cov", "ln_ppp_cov",
                "Hib3_cov", "skilled_ba_cov"))

### covariates for which beta not to shrink ----
(covariates.noshrink = "intercept") #"Iwhoreg6_3")

## fixed hyperparameter to cross-validate ----
rsdlim.fixed = 0.04  # max value of sd_re[j] where N(0, sd_re[j]^2) prior on random effects
lambda.fixed = 0.9  # lambda in laplace for bayesian lasso

## mcmc specifics ----
nMCMC = 1000
nBurn = 1000
adapt_delta_stan = .9



# ######################################################################## #
# quantities required for modeling ----
# ######################################################################## #

## stan object ----
library(rstan)
rstan::rstan_options(auto_write = TRUE)

model.stan_object = rstan::stan_model(file = "VR.stan") #,
                                      #auto_write = rstan::rstan_options("auto_write"=TRUE))



### covariate ----
#head(LASSO_VR_jags$Studies)
#dim(LASSO_VR_jags$Studies)

### misreporting matrix, death counts ----
#head(LASSO_VR_jags$Deaths)
#dim(LASSO_VR_jags$Deaths)

nCause = length(causes) # number of causes
nCovariate = length(covariates) # number of covariates
nCovariate.noshrink = length(covariates.noshrink) # which beta not to shrink

### countries for random effects ----
(countries = sort(unique(studies$isocode)))
nCountry = length(countries)

### studies ----
nStudy = nrow(studies)
studies.name = paste0("study", 1:nStudy)

### design matrix ----
Xmat.shrink = studies[,covariates[!(covariates %in% covariates.noshrink)]]
head(Xmat.shrink)

#### standardize relevant covariates ----
Xmat.shrink = scale(x = Xmat.shrink)
head(Xmat.shrink)

#### appending intercept and not to shrink covariates at the front ----
Xmat = cbind(studies[,covariates.noshrink],
             Xmat.shrink)

# adding intercept to the list
nCovariate = length(covariates)
nCovariate.noshrink = length(covariates.noshrink)

rownames(Xmat) = studies.name
colnames(Xmat) = covariates
head(Xmat)
table(studies$year)


# ######################################################################## #
# model fitting ----
# ######################################################################## #
nchain<- 4

## stan fit ----
ptm0 = Sys.time()
stanfit = rstan::sampling(model.stan_object,
                          data = c(
                            list(
                              "nStudy" = nStudy,
                              "nCause" = nCause,
                              "nCumreport" = c(0,
                                               cumsum(as.numeric(table(deaths$sid)))),
                              "Missreport" = deaths[,causes],
                              "nDeaths" = deaths$n,
                              "nre" = length(countries),
                              "reid" = studies$rG,
                              "K" = nCovariate,
                              "nNoshrink" = nCovariate.noshrink,
                              "Xmat" = Xmat,
                              "sd_betareg_noshrink" = .5,
                              "rsdlim" = rsdlim.fixed,
                              "lambda" = lambda.fixed
                            )
                          ),
                          pars = c(
                            'B', 're', 'sd_re', 'loglik'
                          ),
                          include = T,
                          chains = nchain, 
                          iter = nBurn + nMCMC, warmup = nBurn,
                          control = list('adapt_delta' = adapt_delta_stan),
                          seed = 1)
ptm1 = Sys.time()
ptm1-ptm0 # time taken in seconds


# Save means and SD of predictors for standardisation:
dv <- data.frame(xvar=covariates, mean = colMeans(studies[, covariates], na.rm = T), 
                 sd = apply(studies[, covariates], 2, sd, na.rm=T), 
                 stringsAsFactors = F, row.names = NULL)
## MCMC output ----
MCMCout = rstan::extract(stanfit)

names(MCMCout)


dim(MCMCout$B)
dim(MCMCout$re)
dim(MCMCout$sd_re)

dimnames(MCMCout$B) = list(NULL, covariates, causes[-1])
dimnames(MCMCout$re) = list(NULL,
                            countries,
                            causes[-1])
colnames(MCMCout$sd_re) = causes[-1]

dimnames(MCMCout$B)
dimnames(MCMCout$re)
colnames(MCMCout$sd_re)


## diagnostics ----
### max Rhat ----
(max_Rhat = max(apply(X = MCMCout$B, 2:3,
                      FUN = function(v){
                        
                        rstan::Rhat(v)
                        
                      }),
                apply(X = MCMCout$re, 2:3,
                      FUN = function(v){
                        
                        rstan::Rhat(v)
                        
                      }),
                apply(X = MCMCout$sd_re, 2,
                      FUN = function(v){
                        
                        rstan::Rhat(v)
                        
                      })))

### min bulk ESS ----
(min_ess_bulk = min(apply(X = MCMCout$B, 2:3,
                          FUN = function(v){
                            
                            rstan::ess_bulk(v)
                            
                          }),
                    apply(X = MCMCout$re, 2:3,
                          FUN = function(v){
                            
                            rstan::ess_bulk(v)
                            
                          }),
                    apply(X = MCMCout$sd_re, 2,
                          FUN = function(v){
                            
                            rstan::ess_bulk(v)
                            
                          }))/nMCMC)

### information criterion for model selection performacne ----
#### loo ic ----
MCMCout$loo.out = loo::loo(MCMCout$loglik, 
                           r_eff = loo::relative_eff(exp(MCMCout$loglik), 
                                                     chain_id = rep(1, nrow(MCMCout$loglik))), 
                           cores = 1)

#### waic ----
MCMCout$waic.out = loo::waic(MCMCout$loglik)

ic.df = rbind(MCMCout$waic.out$estimates,
              MCMCout$loo.out$estimates)

ic.df.melt = as.numeric(t(ic.df))
names(ic.df.melt) = paste0(rep(rownames(ic.df), each = ncol(ic.df)),'_',rep(colnames(ic.df), nrow(ic.df)))

MCMCout$ic.df = ic.df.melt

#### mcmc diagnostic ----
(MCMCout$mcmc.diagnostic = c('max_Rhat' = max_Rhat, 'min_ess_bulk' = min_ess_bulk,
                             'num_divergent' = rstan::get_num_divergent(stanfit),
                             'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit)))


## save output ----
if(!dir.exists(output.path)){
  
  dir.create(output.path)
  
}


out = list(MCMCout = MCMCout, 
           time = ptm1-ptm0, rsdlim = rsdlim.fixed,
           lambda = lambda.fixed,  
           nchain = nchain, 
           vdt = causes, 
           vxf = covariates,
           Vars=dv, studies=studies)

save(out, file=paste0("../results/VR_lambda", lambda.fixed, 
                                  "_rsdlim", rsdlim.fixed,"_oldinput.RData"))
#mean(MCMCout$B)

xt<-colMeans(MCMCout$B)
xt

xtr<-colMeans(MCMCout$re)
head(xtr)

load("../results/VR.RData")
xto<-out$output$BUGSoutput$mean$B
xto
head(out$output$BUGSoutput$mean$re)
#out$param$VXF
summary(xt - xto)
out$vxf