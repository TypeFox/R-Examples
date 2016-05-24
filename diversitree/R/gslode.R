make.ode.gslode <- function(info, control) {
  n.var <- info$ny
  n.par <- info$np 
  rtol <- atol <- control$tol
  stepper <- control$gsl.stepper

  ## Some checking:
  if ( length(rtol) != 1 )
    stop("rtol must (currently) be scalar")
  if ( length(atol) != 1 )
    stop("atol must (currently) be scalar")

  time.varying <- isTRUE(info$time.varying)
  if ( time.varying )
    tm <- info$tm

  if ( control$compiled ) {
    model <- info$name.ode
    dll <- info$dll
    derivs <- sprintf("derivs_%s_gslode", model)
    derivs <- getNativeSymbolInfo(derivs, PACKAGE=dll)$address

    if ( time.varying )
      ode <- new(GslOdeTime,     derivs, n.var, tm)
    else
      ode <- new(GslOdeCompiled, derivs, n.var)

  } else {
    derivs <- info$derivs
    ode   <- new(GslOdeR,        derivs, environment(derivs), n.var)
  }

  ## Control parameters (will get tweaked.
  ode$set_control(list(atol=atol, rtol=rtol, algorithm=stepper,
                       hini=1e-4))

  do.set.tm <- time.varying && !control$compiled
  
  function(vars, times, pars) {
    if ( length(pars) != n.par )
      stop("Incorrect parameter length")
    if ( length(vars) != n.var )
      stop("Incorrect variable length")
    if ( length(times) <= 1 )
      stop("Need >= 2 times")
    if ( do.set.tm )
      tm$set(pars)
    ode$run(times, vars, pars)
  }
}
