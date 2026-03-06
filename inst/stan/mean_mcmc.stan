// Analyse a proposed set of prior data, given vague prior pv() represented as a mixture of normals
// Posterior produced from this can be compared with the elicited information
// Quantity of interest is the mean of a continuous outcome.  (For probabilities, see prob_mcmc.stan)

// TESTME

data {
  // Normal mixture components for prior density of ymean
  int<lower=1> nmix;
  vector[nmix] means;
  vector[nmix] sds;
  vector[nmix] probs;

  // Empirical mean and standard error from a sample generated from a normal with mean ymean and some unspecified sampling SD
  // Standard error = sampling SD / sqrt(sample size).  Implicit non-integer sample size OK because we can still "analyse" the "prior data" with this
  real ybar;
  real<lower=0> se;

  int<lower=0,upper=1> pred;
}

parameters {
  real ymean;
  vector[pred] alpha;
}

transformed parameters {
}

model {
  real prior_dens, prior_logdens, prior_logdens_surv;
  prior_dens = 0;
  for (i in 1:nmix){
    prior_dens += exp(log(probs[i]) + normal_lpdf(ymean | means[i], sds[i]));
  }
  prior_logdens = log(prior_dens);
  target += prior_logdens;
  if (pred==0){
    ybar ~ normal(ymean, se);
  } else {
    alpha ~ normal(ymean, 1);
    ybar ~ normal(alpha, se);
  }
}
