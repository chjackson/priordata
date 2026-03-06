## Functions to infer the prior data given elicited info and a sample from pv()
## by fitting a conjugate distribution to pv() 


##' @importFrom betareg betareg
beta_mle_fit <- function(sample_prior, ...){
  bm <- betareg::betareg(sample_prior ~ 1, ...)
  mu <- plogis(coef(bm)[1])
  phi <- exp(coef(bm)[2])
  list(a = mu*phi,
       b = (1-mu)*phi)
}

#' Obtain binomial prior data by conjugacy given elicited quantiles
#' and a sample from the vague prior 
#' 
#' @export
priordata_conjugate_binomial <- function(sample_prior, elicited){
  prior <- beta_mle_fit(sample_prior)
  beta_post <- SHELF::fitdist(vals=elicited$quantiles, probs=elicited$probs, lower=0, upper=1)$Beta
  y <- beta_post["shape1"] - prior$a
  n <- beta_post["shape2"] - prior$b + y
  list(y=round(y), n=round(n))
}


#' @export
normal_mle_fit <- function(sample_prior){ # not actually mle, robust 
  list(mu = median(sample_prior),
       sigma = unname(diff(quantile(sample_prior, c(0.025, 0.975)))) / (qnorm(0.975)-qnorm(0.025)))
}

#' Obtain normal prior data by conjugacy given elicited quantiles
#' and a sample from the vague prior 
#' 
#' @export
priordata_conjugate_normal <- function(sample_prior, elicited){
  prior <- normal_mle_fit(sample_prior)
  if (identical(names(elicited), c("ybar", "se")))
    normal_post <- list(mu=elicited$ybar, sigma=elicited$se)
  else
    normal_post <- SHELF::fitdist(vals=elicited$quantiles, probs=elicited$probs, lower=0, upper=1)$Normal
  npost <- 1/normal_post[["sigma"]]^2
  nprior <- 1 / prior$s^2
  nobs <- npost - nprior
  ybar <- (normal_post[["mu"]]*nobs - nprior*prior$m)/nobs
  se <- sqrt(1/nobs)
  data.frame(ybar=ybar, se=se) # TODO DF OR LIST 
}
