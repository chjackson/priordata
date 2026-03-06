## todo whats the equivs for the non-mixture version
## instead of mixprior, have a stan model and data
## something that we add pdata to.  so a dataset and a Stan model 

#' Posterior from combining binomial prior data with a vague prior
#' approximated by a mixture distribution
#' 
#' @param pdata Binomial prior data: data frame or list with components
#' \code{y} and \code{n}.
#' 
#' @return An \code{rstan} sampling object as returned by \code{sampling()}.
mcmc_mixture_binomial <- function(pdata, mixprior){
  data <- c(list(y=pdata$y, n=pdata$n, 
                 nmix = length(mixprior$probs),
                 pred=0),
            mixprior)
  stanfit <- suppressMessages(
    rstan::sampling(stanmodels$prob_mcmc, data=data,
             refresh=0, open_progress=FALSE, show_messages=FALSE)
  )# TODO split and get diagnostics 
  stanfit
}

#' Posterior from combining normal prior data with a vague prior
#' approximated by a mixture distribution
#'
#' @param pdata Normal prior data: data frame or list with components
#' \code{ybar} and \code{se}.
#' 
#' @return An \code{rstan} sampling object as returned by \code{sampling()}.
mcmc_mixture_normal <- function(pdata, mixprior, ...){
  data <- c(list(ybar=pdata$ybar, se=pdata$se, 
                 nmix = length(mixprior$probs),
                 pred=0),
            mixprior)
  stanfit <- rstan::sampling(stanmodels$mean_mcmc, data=data, ...)
  stanfit
}

mcmc_mixture_extract_estimand <- function(stanfit){
  samdf <- as_draws_df(stanfit)
  estimand <- if ("ymean" %in% names(samdf)) "ymean" else "surv" 
  samdf[[estimand]]
}

#' Squared error loss representing disagreement between elicited
#' information about a proportion and a "prior data" representation
#'
#' Firstly, this function obtains the posterior from combining a vague
#' prior for a logit probability (represented as a mixture of normal
#' distributions) with a set of binomial "prior data" representing an
#' expert's additional knowledge.
#'
#' Secondly, the loss between the posterior and the elicited information
#' is obtained. 
#' 
#' 
#'
loss_mcmc_mixture_binomial <- function(pdata, mixprior, elicited){
  sample_post <- mcmc_mixture_binomial(pdata, mixprior)
  loss_post(sample_post, elicited)
}

loss_mcmc_mixture_normal <- function(pdata, mixprior, elicited){
  sample_post <- mcmc_mixture_normal(pdata, mixprior)
  loss_post(sample_post, elicited)
}

# or could use S3 methods

loss_mcmc_mixture <- function(pdata, mixprior, elicited){
  if (pdata_is_binomial(pdata))
    loss_mcmc_mixture_beta(pdata, mixprior, elicited)
  else if (pdata_is_normal(pdata))
    loss_mcmc_mixture_normal(pdata, mixprior, elicited)
}
