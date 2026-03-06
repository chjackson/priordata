## 2D normal optimisation: ybar and (log) standard error

##' Objective function for searching for optimal prior data that
##' minimises loss between posterior and elicited information
##'
##' Parameterised in `optim` style
##'
##' @param par Prior data as vector transformed to real line:
##'   empirical mean and standard error of implicit normal sample
##'
##' @param mixprior Required if model_fn not specified
##'
##' @param model_fn
##' Function mapping pdata to posterior sample for estimand
##'
##' @export
normal_objective <- function(par, elicited,
                             mixprior=NULL, model_fn=NULL,
                             file = NULL,
                             ...){
  pdata <- list(ybar = par[1],
                se = exp(par[2]))
  sample_post <- mcmc_normal(pdata, mixprior, model_fn, ...)
  res <- loss_post(sample_post, elicited)
  mcse <- mcse_loss(sample_post, elicited)
  pdiv <- attr(sample_post, "pdiv")
  if (!is.null(file)){
    intermediates <- c(pdata$ybar, pdata$se, res, mcse, pdiv)
    cat(paste(intermediates, collapse="  ,  "), "\n", file=file, append=TRUE)
  }
  res
}

mcmc_normal <- function(pdata, mixprior=NULL, model_fn=NULL, ...){
  if (!is.null(model_fn)){
    sample_post <- model_fn(pdata, ...)
  } else {
    stanfit_post <- mcmc_mixture_normal(pdata, mixprior, ...)
    sample_post <- mcmc_mixture_extract_estimand(stanfit_post)
    attr(sample_post, "pdiv") <- get_divergences(stanfit_post)
  }
  sample_post
}

get_divergences <- function(stanfit){
  diags <- rstan::get_sampler_params(stanfit)
  diags <- do.call("rbind", diags)
  pdiv <- mean(diags[,"divergent__"] == 1)
  pdiv
}

normal_opt <- function(pdata_init, elicited,
                       mixprior=NULL, model_fn=NULL,
                       optim_control=NULL, file=NULL, ...){
  par <- c(ybar = pdata_init$ybar,
           logse = log(pdata_init$se))

  opt <- optim(par, fn=normal_objective, elicited=elicited,
               mixprior=mixprior, model_fn=model_fn,
               control=optim_control, file=file, ...)

  par_opt <- opt$par
  pdata_opt <- data.frame(ybar = par_opt[["ybar"]],
                          se = exp(par_opt[["logse"]]))
  pdata_opt
}


normal_opt_bayes <- function(pdata_init, elicited,
                             mixprior=NULL, model_fn=NULL,
                             bc = NULL,
                             file = NULL,
                             ...){
  bc_user <- bc
  bc <- bc_default
  for (i in names(bc_user)) bc[[i]] <- bc_user[[i]]
  obj <- smoof::makeSingleObjectiveFunction(
    name = "normal_objective",
    fn = function(x) {
      normal_objective(par=c(x[["ybar"]],x[["logse"]]),
                       elicited=elicited,
                       mixprior=mixprior, model_fn=model_fn, file=file, ...)
    },
    par.set = ParamHelpers::makeParamSet(
      ParamHelpers::makeNumericParam("ybar", lower=155, upper=165),
      ParamHelpers::makeNumericParam("logse", lower=log(8*0.7), upper=log(8/0.7))
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  des <- ParamHelpers::generateDesign(n = bc$ndes,
                                      par.set = getParamSet(obj),
                                      fun = lhs::randomLHS)
  surr_km <- makeLearner("regr.km", predict.type = "se", covtype = "matern3_2",
                         nugget.estim = TRUE,  control = list(trace = bc$km_trace))
  control <- mlrMBO::makeMBOControl(final.method = "best.predicted",
                                    final.evals = bc$final.evals)
  control <- setMBOControlInfill(control, crit = makeMBOInfillCritAEI())
  if (!is.null(bc$time.budget))
    control <- setMBOControlTermination(control, time.budget = bc$time.budget)
  else if (!is.null(bc$iters))
    control <- setMBOControlTermination(control, iters = bc$iters)
  run <- mbo(obj, design = des, learner = surr_km,
             control = control, show.info = bc$show.info)
  opt <- data.frame(ybar = run$x$ybar,
                    se = exp(run$x$logse),
                    loss = run$y)
  attr(opt, "mbo") <- run
  opt
}

bc_default <- list(
  ndes = 10,
  time.budget = 600,
#  iters = 30,
  final.evals = 10,
  km_trace = FALSE,
  show.info = TRUE
)
