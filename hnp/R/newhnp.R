newhnp <-
function(object, sim=99, conf=.95, halfnormal=T, plot.sim=T, 
         verb.sim=F, how.many.out=F, print.on=F, paint.out=F, 
         col.paint.out, diagfun, simfun, fitfun, ...) {
  
  # checks
  if(missing(fitfun)|missing(simfun)) {
    stop("When running hnp for a non-implemented class or diagnostic measure, you must provide functions for arguments:", "\n", 
         "diagfun (function used for diagnostic extraction - default is resid()),", "\n",
         "simfun (function used for data simulation) and", "\n",
         "fitfun (function used for model fitting).", "\n")
  }
  if(missing(diagfun)) {
    diagfun <- function(obj) residuals(obj)
    cat("Using default residuals computed through residuals().", "\n")
  }
  
  # producing the envelope bands
  n <- length(diagfun(object))
  simulated <- lapply(rep(n, sim), simfun, object)
  if(halfnormal) {
    if(verb.sim) {
      fn <- function(x) {
        out <- sort(abs(diagfun(fitfun(x))))
        cat("Model refitting concluded", "\n")
        return(out)
      }
    } else {
      fn <- function(x) sort(abs(diagfun(fitfun(x))))
    }
  } else {
    if(verb.sim) {
      fn <- function(x) {
        out <- sort(diagfun(fitfun(x)))
        cat("Model refitting concluded", "\n")
        return(out)
      }
    } else {
      fn <- function(x) sort(diagfun(fitfun(x)))
    }
  }
  all.res <- lapply(simulated, fn)
  if(halfnormal) {
    res.original <- sort(abs(diagfun(object)))
  } else {
    res.original <- sort(diagfun(object))
  }
  res <- cbind(res.original, sapply(all.res, cbind))
  
  # now run .makehnp
  .makehnp(obj=res, conf=conf, halfnormal=halfnormal, how.many.out=how.many.out, 
              paint.out=paint.out, col.paint.out=col.paint.out, print.on=print.on, plot.sim=plot.sim, ...)
}
