##' Greedy optimisation procedure to find optimal binomially-distributed
##' prior data representing elicited judgements about a probability
##' 
##' @param pdata_init data frame with components y n giving
##'  initial values
##'
##' @param w half-width of window used to define the neighbourhood of
##' points used in one iteration of the greedy search.  The loss is
##' computed for all points inside this neighbourhood. 
##'
##' @inheritParams ???
##' 
##' @return list with components
##'
##' `loss` data frame with components y, n and loss, giving loss for each y,n pair
##'
##' `sam` data frame of MCMC samples for the parameter of interest, with
##' one row per y, n and MCMC sample.
##'
##' `ecdf` data frame with components y, n, x and ecdf, giving ECDF of the
##'  parameter of interest at the ordinate given by x
##'
##' `opt` (for convenience) data frame with one row, the row of `loss` with
##' the minimum loss
##'
##' TODO switch between mixprior and true
##' consistent interface between binomial and normal 
##' 
##'
##' @md
##' @noRd 
binomial_opt_greedy <- function(pdata_init, model_fn, mixprior, elicited, w=5,
                                loss=TRUE, sam=TRUE, ecdf=TRUE){
  minloss_new <- minloss_old <- Inf 
  pdata_curr <- pdata_init 
  losses_done <- data.frame(y=numeric(), n=numeric(), loss=numeric())
  sam_done <- data.frame(y=numeric(), n=numeric(), sam=numeric())
  ecdf_done <- data.frame(y=numeric(), n=numeric(), x=numeric(), ecdf=numeric())

  while((minloss_new < minloss_old)
        || is.infinite(minloss_new))
  {
    neigh <- neighbourhood(pdata_curr, w)
    neigh_new <- setdiff_df(neigh, losses_done[,c("y","n")])
    res <- pdata_to_sam_binomial(neigh_new, model_fn,
                                        mixprior, elicited)
    losses_done <- rbind(losses_done, res$loss)
    sam_done <- rbind(sam_done, res$sam)
    ecdf_done <- rbind(ecdf_done, res$ecdf)
    opti <- which.min(losses_done$loss)
    minloss_old <- minloss_new
    minloss_new <- losses_done[opti,"loss"]
    pdata_curr <-  losses_done[opti,c("y","n")]
    cat(sprintf("y=%s, n=%s, loss=%s\n", 
                pdata_curr$y, pdata_curr$n, minloss_new))
  }
  list(loss=losses_done,
       sam=sam_done,
       ecdf=ecdf_done,
       opt=losses_done[opti,])
}

##' Like base::setdiff but for data frame rows 
##' @param x, y data frames with same column names 
##' @return rows of x that are not included in y
##' @noRd
setdiff_df <- function(x, y){
  if (is.null(y) || (nrow(y)==0)) return(x)
  xp <- apply(x, 1, paste, collapse=",")
  yp <- apply(y, 1, paste, collapse=",")
  x[!(xp %in% yp),,drop=FALSE]
}

##' Set of 2D pairs of integers within a neighbourhood of a given point
##'
##' @param pdata data frame with one row and components y and n
##' @param window half width 
##' @return data frame with components y and n, and one row
##' for each (y,n) pair with y or n within +- w of the original
##' y and n 
##' @noRd
neighbourhood <- function(pdata, w){
  yrange <- seq(max(pdata$y - w, 0),
                min(pdata$y + w, pdata$n + w))
  nrange <- seq(max(pdata$n - w, 0),
                pdata$n + w)
  dat <- expand.grid(y=yrange, n=nrange)
  dat[dat$y <= dat$n,]
}

#' Calculate loss for a range of proposed binomial data y, n 
#' 
#' @param pdata data frame with one row for each proposed prior dataset
#'
#' @param model_fn Function with one argument, a 1-row data frame or
#'   list with components y and n, and returning a sample from the
#'   posterior given just this dataset and no other data.  This sample
#'   is to be equated with the elicited information.
#'
#' @param mixprior TODO document normal mixture prior specification, to be used
#' if we don't have access to the model. Either model_fn or mixprior should
#' be supplied. 
#'
#' @param elicited TODO
#' 
#' @return tidy data frame with the loss for each proposed dataset given
#' the supplied elicited data
pdata_to_sam_binomial <- function(pdata, model_fn=NULL, mixprior=NULL,
                                  elicited, ecdf_x=NULL, loss="mad"){
  nyn <- nrow(pdata)
  sams <- ecdfs <- lossdf <- vector(nyn, mode="list")
  if (is.null(ecdf_x)) ecdf_x <- ppoints(100)
  for (i in seq_len(nyn)){
    y <- pdata$y[i] 
    n <- pdata$n[i]
    pdi <- list(y=y, n=n)
    sam <- mcmc_binomial(pdi, mixprior, model_fn)
    lossdf[[i]] <- data.frame(y=y, n=n,
                            loss = loss_post(sam, elicited, loss=loss))
    sams[[i]] <- data.frame(y=y, n=n, sam=sam)
    ecdfs[[i]] <- data.frame(y=y, n=n, x=ecdf_x, ecdf = ecdf(sam)(ecdf_x))
  }
  sams <- do.call("rbind", sams)
  ecdfs <- do.call("rbind", ecdfs)
  lossdf <- do.call("rbind", lossdf)
  list(sam=sams, ecdf=ecdfs, loss=lossdf)
}

mcmc_binomial <- function(pdata, mixprior, model_fn){
  if (!is.null(model_fn)){
    sam <- model_fn(pdata)
  } else 
    sam <- mcmc_mixture_binomial(pdata, mixprior)
  sam
}
