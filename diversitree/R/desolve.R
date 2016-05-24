lsoda.trim <- function(...) {
  ret <- t(lsoda(...)[-1,-1,drop=FALSE])
  dimnames(ret) <- NULL
  ret
}

## This sets things up the way that deSolve likes them
derivs.for.deSolve <- function(f)
  function(...) list(f(...))

make.ode.deSolve <- function(info, control) {
  if ( !is.function(info$derivs) )
    stop("info$derivs must be a function")
  derivs <- derivs.for.deSolve(info$derivs)

  rtol <- atol <- control$tol

  if ( isTRUE(info$time.varying) ) {
    tm <- info$tm
    function(vars, times, pars) {
      tm$set(pars)
      lsoda.trim(vars, times, derivs, pars, rtol=rtol, atol=atol)
    }
  } else {
    function(vars, times, pars) {
      lsoda.trim(vars, times, derivs, pars, rtol=rtol, atol=atol)
    }
  }
}
