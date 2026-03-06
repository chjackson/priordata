## These are all unused




#' Conjugate Bayesian inference for a proportion 
#' @param y Outcome count for binomial data
#' @param n Denominator for  binomial data
#' @param aprior Beta prior first shape parameter
#' @param aprior Beta prior second shape parameter
#' @return List with Beta posterior shape parameters
conjugate_binomial <- function(y, n, aprior, bprior){
  list(a = aprior + y,
       b = bprior + n - y)
}

#' Conjugate Bayesian inference for a normal outcome with known variance
#' @param y Observed empirical mean of a sample generated from a normal
#' with unknown mean mu and sampling variance of 1
#' @param se Observed standard error for this empirical mean
#' @param mprior Mean of normal prior for mu
#' @param sprior SD of normal prior for mu
#' @return List with mean and SD of Normal posterior for mu given this prior and data
conjugate_normal <- function(ybar, se, mprior, sprior){
  nprior <- 1 / sprior^2
  nobs <- 1 / se^2
  mpost <- (nprior*mprior + nobs*ybar) / (nprior + nobs)
  spost <- sqrt(1 / (nprior + nobs))
  list(m = mpost, s=spost) # TODO df or list?
}



#' Squared error loss representing disagreement between elicited
#' information about a proportion and a "prior data" representation.
#'
#' Firstly this function obtains the posterior resulting from
#' combining a (typically vague) prior Beta distribution with
#' "prior data" representing an expert's additional knowledge.
#'
#' The sum of squared differences between the posterior probabilities
#' and the elicited probabilities, for a set of quantiles, is then
#' returned.  TODO other losses
#' 
#' @inheritParams conjugate_binomial
#'
#' @param elicited List of elicited information about a proportion, with components
#'
#' `quantiles` Vector of points on [0,1]
#'
#' `probs` Vector of probabilities that the proportion is less than this quantile
#' 
#' @return Sum of squared difference between posterior and elicited probabilities
loss_conj_binomial <- function(pdata, prior, elicited, loss="mad"){
  post <- conjugate_binomial(pdata$y, pdata$n, prior$a, bprior$b)
  lossfn(pbeta(elicited[["quantiles"]], post$a, post$b),
         elicited[["probs"]],
         loss)
}

#' Squared error loss representing disagreement between elicited
#' information about a normal quantity and a "prior data" representation.
#' 
#' @inheritParams conjugate_normal
#' 
#' @param elicited List of elicited information about a normally distributed quantitity, with components
#'
#' `quantiles` Vector of points on the real line
#'
#' `probs` Vector of probabilities that the quantity is less than this quantile
#' 
#' @return Sum of squared difference between posterior and elicited probabilities
loss_conj_normal <- function(pdata, prior, elicited, loss="mad"){
  post <- conjugate_normal(pdata$ybar, pdata$se, prior$m, prior$s)
  lossfn(pnorm(elicited[["quantiles"]], post$m, post$s),
         elicited[["probs"]],
         loss)
}
