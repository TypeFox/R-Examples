do.mle.search.optim <- function(func, x.init, control, lower, upper) {
  control <- modifyList(list(fnscale=-1,
                             ndeps=rep(1e-5, length(x.init)),
                             optim.method="L-BFGS-B"), control)

  optim.method <- control$optim.method
  control.optim <- control[c("fnscale", "ndeps")]

  ans <- optim(x.init, func, method=optim.method,
               control=control.optim, lower=lower, upper=upper)
  names(ans)[names(ans) == "value"] <- "lnLik"  
  ans$optim.method <- optim.method
  
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (optim): ",
            tolower(ans$message))

  ans
}

do.mle.search.subplex <- function(func, x.init, control, lower, upper) {
  ## By default, lower tolerance-- more likely to be met
  control <- modifyList(list(reltol=.Machine$double.eps^0.25,
                             parscale=rep(.1, length(x.init))),
                        control)

  check.bounds(lower, upper, x.init)
  if ( any(is.finite(lower) | is.finite(upper)) )
    func2 <- invert(boxconstrain(func, lower, upper))
  else
    func2 <- invert(func)

  ans <- subplex(x.init, func2, control)
  ans$value <- -ans$value
  names(ans)[names(ans) == "value"] <- "lnLik"

  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (subplex): ",
            tolower(ans$message))
  
  ans
}

do.mle.search.nlminb <- function(func, x.init, control, lower, upper) {
  control.nlminb.ok <- c("eval.max", "iter.max", "trace", "abs.tol",
                         "rel.tol", "x.tol", "step.min")
  control.nlminb <- control[names(control) %in% control.nlminb.ok]
  ans <- nlminb(x.init, invert(func), control=control.nlminb,
                lower=lower, upper=upper)
  names(ans)[names(ans) == "objective"] <- "lnLik"
  names(ans)[names(ans) == "evaluations"] <- "counts"
  ans$lnLik <- -ans$lnLik
  ans <- ans[c("par", "lnLik", "counts", "convergence", "message",
               "iterations")]
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (nlminb): ",
            tolower(ans$message))
  ans
}

do.mle.search.nlm <- function(func, x.init, control, lower, upper) {
  nlm.defaults <-
    list(typsize=rep(1, length(x.init)), print.level=0, ndigit=12,
         gradtol=1e-06, steptol=1e-06, iterlim=100,
         check.analyticals=TRUE)
  control <- modifyList(nlm.defaults, control)
  if ( is.null(control$stepmax) )
    control$stepmax <- max(1000*sqrt(sum((x.init/control$typsize)^2)), 1000)

  ans <- nlm(invert(func), x.init, typsize=control$typsize,
             print.level=control$print.level, ndigit=control$ndigit,
             gradtol=control$gradtol, stepmax=control$stepmax,
             steptol=control$steptol, iterlim=control$iterlim,
             check.analyticals=control$check.analyticals)

  names(ans)[names(ans) == "estimate"] <- "par"
  names(ans)[names(ans) == "minimum"] <- "lnLik"
  names(ans)[names(ans) == "iterations"] <- "counts"
  ans$lnLik <- -ans$lnLik 
  ans <- ans[c("par", "lnLik", "counts", "code", "gradient")]
  
  if ( ans$code > 2 )
    warning("Convergence problems in find.mle: code = ",
            ans$code, " (see ?nlm for details)")  

  ans
}

do.mle.search.minqa <- function(func, x.init, control, lower, upper) {
  if ( !requireNamespace("minqa") )
    stop("This method requires the minqa package")
  
  control <- modifyList(list(minqa.method="newuoa"), control)
  minqa.method <- match.arg(control$minqa.method,
                            c("bobyqa", "newuoa", "uobyqa"))

  control.minqa.ok <- c("npt", "rhobeg", "rhoend", "iprint", "rho",
                        "maxfun")
  control.minqa <- control[names(control) %in% control.minqa.ok]

  check.bounds(lower, upper, x.init)
  if ( any(is.finite(lower) | is.finite(upper)) )
    func2 <- invert(boxconstrain(func, lower, upper))
  else
    func2 <- invert(func)

  opt <- getExportedValue("minqa", minqa.method)
  if ( minqa.method == "bobyqa" )
    ans <- opt(x.init, func2, lower=lower, upper=upper,
               control=control.minqa)
  else
    ans <- opt(x.init, func2, control=control.minqa)
  
  list(par=ans$par, lnLik=ans$fval, counts=ans$feval,
       minqa.method=minqa.method)
}

