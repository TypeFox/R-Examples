lnre <- function (type=c("zm", "fzm", "gigp"),
                  spc=NULL, debug=FALSE,
                  cost=c("chisq", "linear", "smooth.linear", "mse", "exact"),
                  m.max=15, 
                  method=c("Custom", "NLM", "Nelder-Mead", "SANN"),
                  exact=TRUE, sampling=c("Poisson", "multinomial"),
                  bootstrap=0, verbose=TRUE,
                  ...)
{
  type <- match.arg(type)
  sampling <- match.arg(sampling)
  user.param <- list(...)               # model parameters passed by user
  user.pnames <- names(user.param)
  method <- match.arg(method)
  cost <- match.arg(cost)
  
  if (sampling == "multinomial") {
    warning("multinomial sampling not yet implemented, falling back to Poisson sampling")
    sampling <- "Poisson"
  }

  cost.function <- switch(cost,         # implementation of chosen cost function
                          chisq=lnre.cost.chisq,
                          linear=lnre.cost.linear,
                          smooth.linear=lnre.cost.smooth.linear,
                          mse=lnre.cost.mse,
                          exact=lnre.cost.mse) # use MSE cost with adjusted value for m.max

  constructor <- switch(type,           # select appropriate constructor function
                        zm = lnre.zm,
                        fzm = lnre.fzm,
                        gigp = lnre.gigp)

  model <- constructor(param=user.param) # initialize model with user-specifid parameter values
  model$exact <- exact
  model$multinomial <- sampling == "multinomial"
  
  pnames <- names(model$param)          # get list of model parameters

  given.param <- pnames[pnames %in% user.pnames]
  missing.param <- pnames[! (pnames %in% user.pnames)]
  unknown.param <- user.pnames[! (user.pnames %in% pnames)]
  if (length(unknown.param) > 0) warning("unknown parameter(s) ignored: ", unknown.param)

  if (length(missing.param) > 0) {
    ## incomplete model -> estimate parameters from observed frequency spectrum
    if (missing(spc)) stop("parameter(s) ", paste(missing.param, collapse=", ")," not specified")
    if (debug) {
      cat("Estimating parameter(s) <", missing.param, "> from observed spectrum.\n")
      cat("Default values for other parameters:\n")
      print(as.data.frame(model$param))
    }
    ## adjust m.max for "exact" parameter estimation (to match V and first V_m exactly)
    if (cost == "exact") m.max <- max(length(missing.param) - 1, 1)
      
    if (method == "Custom") { # custom estimation uses method call to fall back on default automatically
      model <- estimate.model(model, spc=spc, param.names=missing.param, debug=debug,
                              method=method, cost.function=cost.function, m.max=m.max)
    }
    else {
      model <- estimate.model.lnre(model, spc=spc, param.names=missing.param, debug=debug,
                                   method=method, cost.function=cost.function, m.max=m.max)
    }
    
    model$spc <- spc
    
    if (bootstrap > 0) {
      model$bootstrap <- lnre.bootstrap(model, N(spc), ESTIMATOR=lnre, STATISTIC=identity, replicates=bootstrap, debug=FALSE, verbose=verbose, type=type, cost=cost, m.max=m.max, method=method, exact=exact, sampling=sampling, ...)
    }
  }
  else {
    ## all parameters specified -> no estimation necessary
    if (! missing(spc)) warning("no use for observed frequency spectrum 'spc' (ignored)")
    if (bootstrap > 0) warning("can't bootstrap fully specified model (skipped)")
    ## parameter values have already been set in constructor call above
  }

  model
}
  
