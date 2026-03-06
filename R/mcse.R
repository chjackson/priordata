## Calculate Monte Carlo standard error of an arbitrary summary function
## using nonparametric bootstrap resampling 

mcse_boot <- function(x, summ_fn, B=1000){
  res <- numeric(B)
  n <- length(x)
  for (i in 1:B){
    x_resample <- x[sample(1:n, size=n, replace=TRUE)]
    res[i] <- summ_fn(x_resample)
  }
  sd(res)
}
