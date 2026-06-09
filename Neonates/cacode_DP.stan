
data {
  
  int nStudy; // number of studies
  int<lower=1> K; // total number of covariates including intercept
  int<lower=1,upper=K> nNoshrink; // number of covariates (including intercept) not shrunk by lasso (arranged at the beginning)
  matrix[nStudy,K] Xmat; // covariates. study by covariates
  
    // random effect
  int nre; // number of random effects for each cause
  array[nStudy] int reid; // random effect id for each study
  
    // Vector with first row of deaths for each study
  array[nStudy] int first;
  
   // Vector with last row of deaths for each study
  array[nStudy] int last;

  
  /// DEATHS
  
  int nCause; // number of causes
  int nRows; // number of rows with deaths
  
  // misreporting matrix
  // row: reported causes from studies (\sum_s D_s)
  // column: max number of causes (C)
  // studies are combined along rows
  matrix[nRows,nCause] Missreport;
  
  // number of deaths from reported causes across studies
  // vectorized and concatinated over studies
  array[nRows] int nDeaths;


  // PARAMETERS

  real<lower=0> sd_betareg_noshrink; // N(0, sd_betareg_noshrink^2) prior on beta for unshrunk covariates (for example, intercept)
  real<lower=0> rsdlim; // max value of sd_re[j] where N(0, sd_re[j]^2) prior on random effects
  real<lower=0> lambda; // lambda in laplace for bayesian lasso
  
}

parameters {

  vector<lower=0, upper=rsdlim>[nCause-1] sd_re ; // degree of random effects
  matrix[K, nCause-1] B ; // regression coeffs into covariate by cause matrix
  matrix[nre,nCause-1] re ; // random effects into study by cause matrix

}


model {

  // log-odds matrix as study by cause. WITHOUT reference cause. eq. (5) in mulick et al. (2022)
  matrix[nStudy,nCause-1] XBmat = Xmat*B;

  // loop through studies
  for(s in 1:nStudy){
    // Multinomial model of deaths observed for reported causes
    nDeaths[first[s]:last[s]] ~ multinomial(Missreport[first[s]:last[s],]*softmax(append_row(0, (re[reid[s],] + XBmat[s,])')));

  }

  // priors
  for(j in 1:(nCause-1)){
    re[,j] ~ normal(0, sd_re[j]); // random effect model
    B[1:nNoshrink,j] ~ normal(0, sd_betareg_noshrink); // Normal prior on unshrunk covariates
    B[(nNoshrink+1):K,j] ~ double_exponential(0, 1/lambda); // Laplace(0,𝜆) prior
  }
  sd_re ~ uniform(0, rsdlim); // degree of random effects

}
