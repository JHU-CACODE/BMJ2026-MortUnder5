#######################################################################
# Functions script accompanying 'A Bayesian hierarchical model with integrated 
# covariate selection and misclassification matrices to estimate 
# neonatal and child causes of death' 
#
# by Mulick AR, Oza S, Prieto-Merino D, Villavicencio F, Cousens S, Perin J
#
# BAYESIAN MODEL ESTIMATION
# f.e1:  accepts study (death) input with misclassification matrices and covariates
#
# PREDICT COD DISTRIBUTION (OUT-OF-SAMPLE)
# f.par:  gather model coefficients from MCMC array and other data returned by f.e1
# f.pr2:  Predict COD distributions from covariates (supplied in function) and data returned by f.par
# f.pci2: Calculate credible intervals from data returned by f.pr2
# f.pibs: Calculate PREDICTIONS but doing INDIA manually by state and then put everything together
#
# 1 October 2024
#######################################################################


require(abind)


##################################################Z
####
####   Gather model data from stanfit object, for prediction
####
f.par <- function(MO, NP=200){
  ## MO    Stan object with posterior MCMC coefficient distribution
  ## NP    Number of coefficient sets from which to estimate credible intervals
  # Recover simulation parameters from Stan output
  SA <- rstan::extract(MO$st.output)
  
  # Point estimates of the parameters
  # Means of parameters (add column for referece COD)
  # Betas
  MEB <- cbind(rep(0, dim(SA$B)[2]), apply(SA$B,c(2,3),mean))
  dimnames(MEB) <- list(dimnames(MO$st.data$Xmat)[[2]],
                       dimnames(MO$st.data$Missreport)[[2]])
  # Random effects
  MER <- cbind(rep(0, dim(SA$re)[2]), apply(SA$re,c(2,3),mean))
  dimnames(MER) <- list(MO$st.data$Rname, MO$st.input$vdt)
  
  # Medians of parameters (add column for referece COD)
  # Betas
  Q2B <- cbind(rep(0, dim(SA$B)[2]), apply(SA$B,c(2,3),median))
  dimnames(Q2B) <- list(dimnames(MO$st.data$Xmat)[[2]],
                       dimnames(MO$st.data$Missreport)[[2]])
  # Random effects
  Q2R <- cbind(rep(0, dim(SA$re)[2]), apply(SA$re,c(2,3),median))
  dimnames(Q2R) <- list(MO$st.data$Rname, MO$st.input$vdt)
  
  # Now, NP samples of the parameters
  # selected NP iterations at random
  SS <- sample(dim(SA$B)[1], NP) 
  # prepare array of fixed effects
  BM <- SA$B[SS,,]
  # add 0 coefficients for reference cause of death
  BM <- abind(array(0, dim=dim(BM)[1:2]), BM)
  # Add names to matrix
  dimnames(BM) <- list(c(1:NP),
                       dimnames(MO$st.data$Xmat)[[2]],
                       dimnames(MO$st.data$Missreport)[[2]])
  # prepare array of random effects
  RM <- SA$re[SS,,]
  # Add 0 random effect for reference cause of death
  RM <- abind(array(0, dim=dim(RM)[1:2]), RM)
  # Add names to matrix
  dimnames(RM) <- list(c(1:NP),
                       MO$st.data$Rnames,
                       dimnames(MO$st.data$Missreport)[[2]])
  # prepare array of RESD
  SM <- SA$sd_re[SS,]
  # Add 0 random effect SD for reference cause of death
  SM <- cbind(rep(0, NP), SM)
  # Add names to matrix
  dimnames(SM) <- list(c(1:NP),
                       dimnames(MO$st.data$Missreport)[[2]])
  return(list(BM=BM, RM=RM, SM=SM, MEB=MEB, MER=MER, Q2B=Q2B, Q2R=Q2R, st.input=MO$st.input, st.data=MO$st.data))
}

###
###  END OF FUNCTION
###
################################################ #



################################################ #
####
####   Predictions with an array of coefficients from f.par()
####   
f.pr2 <- function(PA, NPD, ID="pid", PE="any", PLW=1){
  
  ## PA   object produced with function f.par()
  ## NPD  Name of data set with covariates to be used in the prediction
  ## ID   variable in NPD that uniquely identify observations
  ## PE   Period label ("early","late","any")
  ## PLW  distinguishes between preterm and low birth weight as COD (1/0)
  
  ## This function makes predictions with fixed effects only and then adds 
  ## ONE random effect term selected at random among some candidates in RT in each MCMC iteration.
  # Prepare prediction data
  VXN <- names(PA$st.data$xmeans)  # vector of numerical covariates
  VXF <- dimnames(PA$BM)[[2]]  # vector of all covariates
  S   <- dim(get(NPD))[1]  # number of data points to predict
  K   <- length(VXF) # number of covariates including intercept
  H   <- length(VXN) # number of numerical covariates scaled
  C   <- dim(PA$BM)[3] # number of causes of death
  N   <- dim(PA$BM)[1] # number of simulations
  RISO <- dimnames(PA$RM)[[2]] # Countries with random effects in the model
  # Prepare raw variables dataset from prediction sample
  DX <-mutate(get(NPD), per.early=ifelse(PE=="early",1,0),
              per.late=ifelse(PE=="late",1,0),
              premvslbw=PLW) %>%
    dplyr::select(all_of(c(ID, "iso3", "country", "year", VXF)))
  # Scale the numerical columns with means and SD from model
  DX[,VXN] <- scale(DX[,VXN], PA$st.data$xmeans, PA$st.data$xsd)
  
  # Point estimates with the MEANS of the beta-coefficients (fixed effects)
  POE <- cbind(DX[,1:4], as.matrix(DX[,VXF]) %*% PA$MEB) %>% 
    pivot_longer(cols=colnames(PA$MEB), names_to="cod", values_to="lf.me")
  # Add means of random effects  
  POE <- as.data.frame(PA$MER) %>% rownames_to_column(var="iso3") %>%
    pivot_longer(cols=colnames(PA$MEB), names_to="cod", values_to="lr.me") %>%
    left_join(POE,.)
  # Add Point estimates with the MEDIANS of the coefficients (fixed effects)
  POE <- cbind(DX[,1:4], as.matrix(DX[,VXF]) %*% PA$Q2B) %>% 
    pivot_longer(cols=colnames(PA$Q2B), names_to="cod", values_to="lf.q2") %>%
    left_join(POE,.)
  POE <- as.data.frame(PA$Q2R) %>% rownames_to_column(var="iso3") %>%
    pivot_longer(cols=colnames(PA$Q2B), names_to="cod", values_to="lr.q2") %>%
    left_join(POE,.) %>%
    mutate(lr.me=ifelse(is.na(lr.me),0,lr.me),
           lr.q2=ifelse(is.na(lr.q2),0,lr.q2),
           lfr.me=(lf.me+lr.me), lfr.q2=(lf.q2+lr.q2)) %>%
    group_by(pid, iso3, country, year) %>%
      mutate(pf.me=exp(lf.me)/sum(exp(lf.me)),
             pf.q2=exp(lf.q2)/sum(exp(lf.q2)),
             pr.me=exp(lfr.me)/sum(exp(lfr.me)),
             pr.q2=exp(lfr.q2)/sum(exp(lfr.q2))) %>%
    ungroup()

  REX  <- unique(DX$iso3) # countries that do not have random effects
  REXi <- REX[REX %in% RISO] # countries within estimation sample
  REXo <- REX[!(REX %in% RISO)] # countries NOT in estimation sample
  # loop through number of mcmc samples
  for (i in 1:N){
    # First calculate random effects of countries in estimation sample
    LOi <- rownames_to_column(as.data.frame(PA$RM[i,REXi,]), var="iso3") %>%
      # add randomly drawn random effects for countries not in estimation sample
      bind_rows(cbind(data.frame(iso3=REXo), mvrnorm(length(REXo), rep(0,C), cov(PA$RM[i,,])))) %>% 
      pivot_longer(cols=dimnames(PA$RM)[[3]], names_to="cod", values_to="ref")
    # Calculate fixed effects and add the random effects:
    LOi <- cbind(DX[,1:4], as.matrix(DX[,VXF]) %*% PA$BM[i,,]) %>% 
      pivot_longer(cols=dimnames(PA$BM)[[3]], names_to="cod", values_to="fef") %>%
      left_join(LOi) %>% mutate(sample=i)
    if(i==1) LO <- LOi else LO <- bind_rows(LO, LOi) 
  }
  LO <- mutate(LO, fref=fef+ref) %>%
    group_by(pid, iso3, country, year, sample) %>% 
    mutate(pf=exp(fef)/sum(exp(fef)), pr=exp(fref)/sum(exp(fref))) %>%
    ungroup()
  # return:
  return(list(Point_estimates=POE, Predictions=LO[,c("pid","iso3","country", "year","cod", "sample","pf","pr")], Prediction.Data=get(NPD), st.input=PA$st.input))
}


###
###  END OF FUNCTION
###
################################################ #






## Previous versions of functions




################################################ #
####
####   Predictions with an array of coefficients from f.par()
####   
# f.pr2 <- function(PA, NPD, ID="pid", PE="any", PLW=1){
#   
#   ## PA   object produced with function f.par()
#   ## NPD  Name of data set with covariates to be used in the prediction
#   ## ID   variable in NPD that uniquely identify observations
#   ## PE   Period label ("early","late","any")
#   ## PLW  distinguishes between preterm and low birth weight as COD (1/0)
#   
#   ## This function makes predictions with fixed effects only and then adds 
#   ## ONE random effect term selected at random among some candidates in RT in each MCMC iteration.
#   # Prepare prediction data
#   VXN <- names(PA$st.data$xmeans)  # vector of numerical covariates
#   VXF <- dimnames(PA$BM)[[2]]  # vector of all covariates
#   S   <- dim(get(NPD))[1]  # number of data points to predict
#   K   <- length(VXF) # number of covariates including intercept
#   H   <- length(VXN) # number of numerical covariates scaled
#   C   <- dim(PA$BM)[3] # number of causes of death
#   N   <- dim(PA$BM)[1] # number of simulations
#   RISO <- dimnames(PA$RM)[[2]] # Countries with random effects in the model
#   
#   # Prepare raw variables dataset from prediction sample
#   DX <-mutate(get(NPD), per.early=ifelse(PE=="early",1,0),
#               per.late=ifelse(PE=="late",1,0),
#               premvslbw=PLW) %>%
#     dplyr::select(all_of(c(ID, "iso3","year", VXF)))
#   # Scale the numerical columns with means and SD from model
#   DX[,VXN] <- scale(DX[,VXN], PA$st.data$xmeans, PA$st.data$xsd)
#   
#   # Extend Random Effects Matrix to countries in prediction data (RX)
#   REX <- unique(DX$iso3)
#   REX <- REX[!(REX %in% dimnames(PA$RM)[[2]])]
#   RX <- array(rnorm(length(rep(c(PA$SM),length(REX))), 0, rep(c(PA$SM),length(REX))), 
#               dim=c(dim(PA$SM), length(REX)),
#               dimnames=list(dimnames(PA$SM)[[1]], dimnames(PA$SM)[[2]], REX))
#   RX <- aperm(RX, c(1,3,2))
#   # array combining random effect for in-sample and out-of-sample iso3
#   RT <- abind(PA$RM, RX, along = 2)
#   
#   # # Prepare matrix to store predictions
#   LF <- array(NA, dim = c(N,S,C))  # matrix of fixed Predictions
#   dimnames(LF) <- list(dimnames(PA$BM)[[1]], get(NPD)[[ID]], dimnames(PA$BM)[[3]])
#   LR <- LF  # matrix of random Predictions
#   for (i in 1:N){
#     LF[i,,] <- as.matrix(DX[,VXF]) %*% PA$BM[i,,] # logodds with fixed effects
#     LR[i,,] <-  LF[i,,] + RT[i,match(get(NPD)[["iso3"]], dimnames(RT)[[2]]),]
#   }  
#   PF <- aperm(apply(LF, c(1,2), function(x) exp(x)/sum(exp(x))), perm=c(2,3,1))
#   PR <- aperm(apply(LR, c(1,2), function(x) exp(x)/sum(exp(x))), perm=c(2,3,1))
#   # return:
#   return(list(PF=PF, PR=PR, Prediction.Data=get(NPD), st.input=PA$st.input))
# }
# 

###
###  END OF FUNCTION
###
################################################ #



############################################## #
###
###   CALCULATE CIs from PREDICTIONS returned by f.pr2
###
# f.pci2 <- function(PM, PR=c(2.5, 25, 50, 75, 97.5)){
#   ## PM   object from function f.pr2()
#   ## PR   Posterior percentiles of predicted fractions: c(2.5, 25, 50, 75, 97.5)
#   PP <- NULL
#   for(c in dimnames(PM$PF)[[3]]){
#     xf <- data.frame(t(apply(PM$PF[,,c], 2, function(x) c(p_me=mean(x, na.rm=T), p_sd=sd(x, na.rm=T), quantile(x, PR/100, na.rm=T))))) %>%
#       rownames_to_column(var="pid") %>% mutate(cause=c, type="fixed")
#     xr <- data.frame(t(apply(PM$PR[,,c], 2, function(x) c(p_me=mean(x, na.rm=T), p_sd=sd(x, na.rm=T), quantile(x, PR/100, na.rm=T))))) %>%
#       rownames_to_column(var="pid") %>% mutate(cause=c, type="random")
#     if(is.null(PP)) PP <- bind_rows(xf,xr) else PP <- bind_rows(PP,xf,xr)
#   }
#   names(PP)[4:(length(PR)+3)] <- paste0("p_",PR)
#   PP <- left_join(PM$Prediction.Data, mutate(PP, pid=as.numeric(pid)), by="pid")
#   return(list(Pred=PP, st.input=PM$st.input))
# }
###
###  END OF FUNCTION
###
################################################ #







