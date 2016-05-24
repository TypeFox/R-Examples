ucminf = function(par, fn, gr = NULL, ..., control = list(), hessian=0) {
  con <- list(trace=0, grtol=1e-6, xtol=1e-12, stepmax=1, maxeval=500,
              grad='forward',gradstep=c(1e-6,1e-8), invhessian.lt = NULL)
  stopifnot(names(control) %in% names(con))
  con[(namc <- names(control))] <- control
  stopifnot(length(con$gradstep)==2,con$grad %in% c('forward','central'))
  fnstr <- quote(fn(.x, ...))
  grstr <- quote(gr(.x, ...))
  rho = new.env(parent = environment())
  n <- length(par)
  eps <- c(con$grtol,con$xtol)
  if(!is.null(gr)) { grad <- 0 }
  else { grad <- ifelse(con$grad=='forward',1,2) } #else central
  iw <- n*ceiling(max(n+1,(n+11)/2)) + 10
  w <- rep(0,iw)
  trace <- con$trace>0
  icontr = 1+trace+2*!is.null(con$invhessian.lt)
  if(!is.null(con$invhessian.lt))
    w[(4*n+1):(4*n+n*(n+1)/2)] <-  con$invhessian.lt #con$invhessian[logicMat]
  par0 <- rep(0,n)
  for(i in 1:n) #avoid that par from calling env. is overwritten
    par0[i] = par[i]
  xname <- as.double(rep(0,n))
  names(xname) <- names(par)
  assign(".f",      fnstr                   , envir = rho)
  assign(".gr",     grstr                   , envir = rho)
  assign(".n",      as.integer(n)           , envir = rho)
  assign(".x",      xname                   , envir = rho)
  assign(".par",    as.double(par0)         , envir = rho)
  assign(".stepmax",as.double(con$stepmax)  , envir = rho)
  assign(".eps",    as.double(eps)          , envir = rho)
  assign(".maxfun", as.integer(con$maxeval) , envir = rho)
  assign(".w",      as.double(w)            , envir = rho)
  assign(".iw",     as.integer(iw)          , envir = rho)
  assign(".icontr", as.integer(icontr)      , envir = rho)
  assign(".grad",   as.integer(grad)        , envir = rho)
  assign(".grstep", as.double(con$gradstep) , envir = rho)
  #
  .Call(mfopt, rho)
  #
  W <- get(".w", envir = rho)
  icontr <- get(".icontr", envir = rho)
  ans = list(
    par = get(".par", envir = rho),
    value = W[1],
    convergence = icontr,
    message = switch(as.character(icontr),
      '1' ='Stopped by small gradient (grtol).',
      '2' ='Stopped by small step (xtol).',
      '3' ='Stopped by function evaluation limit (maxeval)',
      '4' ='Stopped by zero step from line search',
      '-2'="Computation did not start: length(par) = 0.",
      '-4'="Computation did not start: stepmax is too small.",
      '-5'="Computation did not start: grtol or xtol <= 0.",
      '-6'="Computation did not start: maxeval <= 0.",
      '-7'="Computation did not start: given hessian not pos. definite.",
      '-8'="Computation did not start: iw too small."
      )
    )
  if(0<icontr) {
    if(hessian == 1) {
      if(suppressPackageStartupMessages(suppressWarnings(require("numDeriv")))) {
        p0 <- ans$par
        names(p0) <- names(par)
        ans$hessian <- hessian(fn, p0, method = "Richardson", ...)
      } else {
        cat("Skipped hessian estimation - package 'numDeriv' must be installed for hessian option 1", sep = "\n")
      }
    }
    if(hessian == 2 | hessian == 3) {
      logicMat <- (matrix(-(1:n^2),n,n,byrow=TRUE)+matrix(1:n^2,n,n))<=0
      COV <- matrix(0,n,n)
      COV[logicMat] <- W[(4*n+1):(4*n+n*(n+1)/2)]
      COV <- t(COV)+COV-diag(diag(COV))
      ans$invhessian <- COV
    }
    if(hessian == 3)
      ans$hessian <- solve(COV)
    ans$invhessian.lt <- W[(4*n+1):(4*n+n*(n+1)/2)]
    ans$info = c( maxgradient = W[2],
                  laststep    = W[3],
                  stepmax     = get(".stepmax", envir = rho),
                  neval       = get(".maxfun", envir = rho)
                )
  }
  if(trace) {
    cat(paste(ans$message,'\n'))
    if(!is.null(ans$info))
      print(ans$info)
  }
  nm <- names(par)
  if (!is.null(nm)) {
    names(ans$par) <- nm
    if(!is.null(ans$hessian))
      colnames(ans$hessian) <- rownames(ans$hessian) <- nm
  }
  ans
}
