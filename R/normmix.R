
dnormmix <- function(x, p, mu, sigma){
  res <- rep(0, length(x))
  for (i in seq_along(p)){
    res <- res + p[i]*dnorm(x, mu[i], sigma[i])
  }
  res
}

dnormmix_mp <- function(x, mp){
  res <- rep(0, length(x))
  for (i in seq_along(mp$probs)){
    res <- res + mp$probs[i]*dnorm(x, mp$mean[i], mp$sds[i])
  }
  res
}

##' @importFrom flexmix prior parameters
dnormmix_fl <- function(x, fl){
  p <- flexmix::prior(fl)
  mu <- flexmix::parameters(fl)[1,]
  sigma <- flexmix::parameters(fl)[2,]
  dnormmix(x, p, mu, sigma)
}


##' @noRd
##' @export
fl_mixprior <- function(fl){
  list(probs = as.array(prior(fl)),
       means = as.array(parameters(fl)[1,]),
       sds = as.array(parameters(fl)[2,]))
}
