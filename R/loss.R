
#' Loss representing disagreement between elicited information about
#' an arbitrary quantity and a sample from a posterior resulting from
#' a "prior data" representation.
#'
#' @param sample_post Sample from the posterior resulting from combining
#' a vague prior with "prior data" representing an expert's additional knowledge
#'
#' @param elicited List of elicited information about a quantity, with components
#'
#' `quantiles` Vector of possible values for the quantity
#'
#' `probs` Vector of probabilities that the quantity is less than this quantile
#'
#' @return Sum of squared difference between posterior and elicited probabilities
#' @noRd 
loss_post <- function(sample_post, elicited, loss="mad"){
  cdf_sample <- ecdf(sample_post)(elicited$quantiles)
  cdf_elic <- elicited$probs
  lossfn(cdf_sample, cdf_elic, loss)
}

mcse_loss <- function(sample_post, elicited, loss="mad", B=1000){
  fn <- function(x) loss_post(x, elicited, loss="mad")
  mcse_boot(sample_post, fn, B=B)
}

loss_stanfit <- function(stanfit, elicited, estimand=NULL, loss="mad"){
  sample_post <- mcmc_extract_estimand(stanfit, estimand)
  data.frame(
    loss = loss_post(sample_post, elicited, loss),
    mcse = mcse_loss(sample_post, elicited, loss)
  )
}

mcmc_extract_estimand <- function(stanfit, estimand=NULL){
  samdf <- posterior::as_draws_df(stanfit)
  if (is.null(estimand))
    estimand <- if ("ymean" %in% names(samdf)) "ymean" else "surv"
  samdf[[estimand]]
}

lossfn <- function(proposed, true, loss){
  switch(loss,
         "mse" = mean((proposed - true)^2),
         "rmse" = sqrt(mean((proposed - true)^2)),
         "mad" = mean(abs(proposed - true))
         )
}
