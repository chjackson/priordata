default_priors <- c(ms0=0, ss0=0.5,
                    mds=0, sds=0.5,
                    mh1=log(100), sh1=0.7,
                    ph0a=1, ph0b=1,
                    mingamma=0, maxgamma=20,
                    msd=0, ssd=0.2) # assume height measured accurately. andrew used.

basepars <- c("h0","h1","s0","s1","gamma","sd_height")

no_data <- list(n=0,
                height=as.array(numeric()),
                age=as.array(numeric()))

no_elic <- list(n_elic=0, age_elic=as.array(numeric()),
                ybar_elic=as.array(numeric()), se_elic=as.array(numeric()),
                pred=0)

height_prior_sample <- function(n=1, priors=NULL, age=NULL){
  p <- as.list(default_priors)
  for (i in seq_along(priors)) p[[i]] <- priors[[i]]
  h1 <- exp(rnorm(n, p$mh1, p$sh1))
  ph0 <- exp(rbeta(n, p$ph0a, p$ph0b))
  h0 <- h1*ph0
  s0 <- exp(rnorm(n, p$ms0, p$ss0))
  ds <- exp(rnorm(n, p$mds, p$sds))
  s1 <- s0 + ds
  gamma <- runif(n, p$mingamma, p$maxgamma)
  if (is.null(age)) age <- seq(0, 30, by=0.1)
  PB1(age, h0, h1, s0, s1, gamma)
}

#' @export
PB1 <- function(age, h0, h1, s0, s1, gamma){
  h1 - 2*(h1 - h0) / (exp(s0*(age - gamma)) + exp(s1*(age - gamma)))
}

#' @param obsdata list or data frame with components age and height
#' @noRd
#' @export
height_standat <- function(obsdata,
                           priors=NULL,
                           elicdata=NULL,
                           age_pred = seq(0, 25, by=0.5),
                           pred = TRUE
                           ){
  if (is.null(obsdata)) {
    standat <- no_data
  } else {
#  standat$age <- (standat$age - mean(standat$age)) / sd(standat$age)
#  standat$height <- (standat$height - mean(standat$height)) / sd(standat$height)
    standat <- as.list(obsdata)
    standat$min_gamma <- min(standat$age)
    standat$max_gamma <- max(standat$age)
    standat$n <- length(standat$age)
  }

  priors_use <- default_priors
  for (i in seq_along(priors))
    priors_use[names(priors)[[i]]] <- priors[[i]]
  standat <- c(standat, priors_use)

  standat$n_pred <- length(age_pred)
  standat$age_pred <- as.array(age_pred)

  if (is.null(elicdata))
    elicdata_use <- no_elic
  else {
    elicdata_use <- list(
      age_elic = as.array(elicdata$age),
      ybar_elic = as.array(elicdata$ybar),
      se_elic = as.array(elicdata$se),
      n_elic = length(elicdata$age),
      pred = as.numeric(pred)
    )
  }
  standat <- c(standat, elicdata_use)

  standat
}


#' Function that maps prior data, and age that it refers to, to a sample
#' from the posterior for mean age.  Priors pv from default_priors
height_model_fn <- function(pdata, agei, ...){
  pdata$age <- agei
  standat_y <- height_standat(obsdata=NULL, elicdata=pdata, age_pred = agei)
  stanfit <- sampling(stanmodels$height, data=standat_y, ...)

  sample_post <- posterior::extract_variable(stanfit, "height_pred[1]")
  attr(sample_post, "diags") <- c(
    age = agei,
    rhat = posterior::rhat(sample_post),
    pdiv = get_divergences(stanfit)
  )
  sample_post
}
