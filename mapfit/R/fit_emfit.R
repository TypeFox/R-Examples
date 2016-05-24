#### emfit

emfit.options <- function() {
  list(maxiter = 2000,
    reltol = sqrt(.Machine$double.eps),
    abstol = +Inf)
}

emfit.verbose <- function() {
  list(emstep = FALSE,
    emprogress = 1)
}


emfit <- function(model, data, initialize = TRUE, control = list(), verbose = list(), ...) {
  con <- emfit.options()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control

  ver <- emfit.verbose()
  nmsC <- names(ver)
  ver[(namc <- names(verbose))] <- verbose

  conv <- FALSE

  if (initialize)
  model <- emfit.init(model, data, verbose, ...)

  eres0 <- emfit.estep(model, data, ...)
  model0 <- emfit.mstep(model, eres0$eres, data, ...)
  for (iter in 1:con$maxiter) {
    eres1 <- emfit.estep(model0, data, ...)

    if (!is.finite(eres1$llf)) {
      warning(gettextf("LLF becomes Inf/NaN/NA at %d", iter))
      error <- c(Inf, Inf)
      model1 <- model0
      break
    }

    model1 <- emfit.mstep(model0, eres1$eres, data, ...)

    error <- c(abs(eres1$llf - eres0$llf), abs((eres1$llf - eres0$llf)/eres0$llf))

    if (eres1$llf - eres0$llf < 0)
    warning(gettextf("LLF decreses: iter=%d llf.diff=%e", iter, eres0$llf - eres1$llf))

    if (ver$emstep)
      if (iter %% ver$emprogress == 0)
    message(gettextf(" iter=%d llf=%e (aerror=%e, rerror=%e)", iter, eres1$llf, error[1], error[2]))

    if ((error[1] < con$abstol) && (error[2] < con$reltol)) {
      conv <- TRUE
      break
    }

    eres0 <- eres1
    model0 <- model1
  }

  list(model=model1,
    llf=eres1$llf,
    iter=iter,
    convergence=conv,
    aerror=error[1],
    rerror=error[2],
    df=emfit.df(model1, ...),
    aic=-2*(eres1$llf-emfit.df(model1, ...)),
    control=con,
    verbose=ver)
}

