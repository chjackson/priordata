data {
  // Normal mixture components for prior density of log(-log(survival)) = log(cumulative hazard)
  int<lower=1> nmix;
  vector[nmix] means;
  vector[nmix] sds;
  vector[nmix] probs;

  // Count and denominator
  int<lower=0> y;
  int<lower=0> n;
}

parameters {
  real logH; // log cumulative hazard
}

transformed parameters {
  real surv;
  surv = exp(-exp(logH)); // survival probability
}

model {
  real prior_dens, prior_logdens, prior_logdens_surv;
  prior_dens = 0;
  for (i in 1:nmix){
    prior_dens += exp(log(probs[i]) + normal_lpdf(logH | means[i], sds[i]));
  }
  prior_logdens = log(prior_dens);
  target += prior_logdens;
  y ~ binomial(n, surv);
}
