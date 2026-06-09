
data {

  int nStudy; // number of studies
  int nCause; // number of causes
  
  // cumulative frequencies of reported causes from studies
  // starts at 0
  // helps us track study-specific death counts and causes
  array[nStudy + 1] int nCumreport;
  
  // misreporting matrix
  // row: reported causes from studies (\sum_s D_s)
  // column: max number of causes (C)
  // studies are combined along rows
  matrix[nCumreport[nStudy+1],nCause] Missreport;
  
  // number of deaths from reported causes across studies
  // vectorized and concatinated over studies
  array[nCumreport[nStudy+1]] int nDeaths;
  
  // random effect
  int nre; // number of random effects for each cause
  array[nStudy] int reid; // random effect id for each study

  int<lower=1> K; // total number of covariates including intercept
  int<lower=1,upper=K> nNoshrink; // number of covariates including intercept (arranged at the beginning)
  matrix[nStudy,K] Xmat; // covariates. study by covariates

  real<lower=0> sd_betareg_noshrink; // N(0, sd_betareg_noshrink^2) prior on beta for unshrunk covariates (for example, intercept)
  real<lower=0> rsdlim; // max value of sd_re[j] where N(0, sd_re[j]^2) prior on random effects
  real<lower=0> lambda; // lambda in laplace for bayesian lasso
  
}

parameters {

  vector[K*(nCause-1)] betareg_v; // vectorized (looping over covariate) regression coeffs
  vector[nre*(nCause-1)] re_v; // vectorized (looping over cause) random effects
  vector<lower=0, upper=rsdlim>[nCause-1] sd_re; // degree of random effects

}

transformed parameters {

  matrix[K,nCause-1] B = to_matrix(betareg_v, nCause-1, K)'; // regression coeffs into covariate by cause matrix
  matrix[nre,nCause-1] re = to_matrix(re_v, nre, nCause-1); // random effects into study by cause matrix
  
  // log-odds matrix as study by cause. WITHOUT reference cause. eq. (5) in mulick et al. (2022)
  matrix[nStudy,nCause-1] XBmat = Xmat*B;
  
  vector[nStudy] loglik; // log-likelihood
  
  // loop through studies
  for(s in 1:nStudy){
    
    // Multinomial model of deaths observed for reported causes
    loglik[s] = multinomial_lpmf(nDeaths[(nCumreport[s]+1):nCumreport[s+1]] | Missreport[(nCumreport[s]+1):nCumreport[s+1],]*softmax(append_row(0, (re[reid[s],] + XBmat[s,])')));

  }

}

model {

  // prior
  // random effect model
  for(j in 1:(nCause-1)){

    target += normal_lpdf(re_v[((j-1)*nre + 1):(j*nre)]/sd_re[j] | 0, 1) - log(sd_re[j]); // for each cause

  }
  target += uniform_lpdf(sd_re | 0, rsdlim); // degree of random effects

  // prior on beta
  // N(0, sd_betareg_noshrink^2) prior on unshrunk covariates (for example, intercept)
  target += normal_lpdf(betareg_v[1:(nNoshrink*(nCause-1))] | 0, sd_betareg_noshrink);

  // Laplace(0,𝜆) prior on unshrunk covariates (everything except that are shrunk)
  target += double_exponential_lpdf(betareg_v[(nNoshrink*(nCause-1) + 1):(K*(nCause-1))] | 0, 1/lambda);

  // log-likelihood
  target += loglik;

}
