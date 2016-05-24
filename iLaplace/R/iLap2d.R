##' @name iLap2d
##' @title Improved Laplace approximation for bivariate integrals of unimodal functions
##' @description This function is similar to \code{iLap} except that it handles only bivariate integrals of user-written unimodal functions.
##' @usage iLap2d(fullOpt, ff, ff.gr, ff.hess,
##'          control = list(sp.points = 100, delta = 15), ...)
##'
##' @param fullOpt A list containing the minium (to be accesed via \code{fullOpt$par}), the value of the function at the minimum (to be accessed via \code{fullOpt$objective}) and the Hessian matrix at the minimum (to be accessed via \code{fullOpt$hessian}
##' @param ff The minus logarithm of the integrand function (the \code{h} function, see Details).
##' @param ff.gr The gradient of \code{ff}, having the exact same arguments as  \code{ff}
##' @param ff.hess The Hessian matrix of\code{ff}, having the exact same arguments as  \code{ff}
##' @param control A named list of control parameters with elements \code{sp.points}, \code{delta} and \code{n.cores}. \code{sp.points} sets the number points for the spline evaluations; \code{delta} controls the length of the inteval of integration.
##' @param ... Additional arguments to be passed to \code{ff}, \code{ff.gr} and \code{ff.hess}
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}
##'
##' @examples
##' # The negative integrand function in log
##' # is the negative log-density of the multivariate
##' # Student-t density centred at 0 with unit scale matrix
##' ff <- function(x, df) {
##'        d <- length(x)
##'        S <- diag(1, d, d)
##'        S.inv <- solve(S)
##'        Q <- colSums((S.inv %*% x) * x)
##'        logDet <- determinant(S)$modulus
##'        logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) +
##'        logDet) - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
##'        return(-logPDF)
##'        }
##'
##' # the gradient of ff
##' ff.gr <- function(x, df){
##'             m <- length(x)
##'             kr = 1 + crossprod(x,x)/df
##'             return((m+df)*x/(df*kr))
##'             }
##'
##' # the Hessian matrix of ff
##' ff.hess <- function(x, df) {
##' m <- length(x)
##' kr <- as.double(1 + crossprod(x,x)/df)
##' ll <- -(df+m)*2*tcrossprod(x,x)/(df*kr)^2.0
##' dd = (df+m)*(kr - 2*x^2/df)/(df*kr^2.0)
##' diag(ll) = dd;
##' return(ll)
##' }
##'
##' df = 5
##' opt <- nlminb(rep(1,2), ff, gradient = ff.gr, hessian = ff.hess, df = df)
##' opt$hessian <- ff.hess(opt$par, df = df);
##'
##' # improved and standard Laplace approximation
##' # true value is 0.0
##' iLap <- iLap2d(opt, ff = ff, ff.gr = ff.gr, ff.hess = ff.hess,
##'              df = df)
##' Lap <- log(2*pi) - opt$objective - 0.5*determinant(opt$hessian)$mod;
##'
##' @return a double, the logarithm of the integral
##' @export
iLap2d <- function(fullOpt, ff, ff.gr, ff.hess,
                     control = list(sp.points = 100, delta = 15), ...)
{
  if(is.null(control$sp.points)) control$sp.points <- 100
  if(is.null(control$delta)) control$delta <- 15

  m = length(fullOpt$par)
  se = SEv(fullOpt$hessian, m)
  oo = seqMat(par = fullOpt$par, se, lengthOut = control$sp.points,
              q = m, delta = control$delta)
  par.val <- oo$parVal

  fullOpt$ldblock =  ldetHessBlocks(fullOpt$hessian, m)


  log.ncost = log.den = 0.0
  lo = oo$lo
  up = oo$up
  par.val <- par.val[-(control$sp.points + 1),]

  # the marginal
  marg2 = function(x1, ...) {
    out <- tryCatch({
      tmp = function(x) ff(c(x1, x), ...)
      gr.tmp = function(x) ff.gr(c(x1, x), ...)[-1]
      hes.tmp = function(x) ff.hess(c(x1,x), ...)[-1,-1]
      tmpOpt = nlminb(fullOpt$par[2], tmp, gradient = gr.tmp)
      tmpOpt$hessian = hes.tmp(tmpOpt$par)
      ans = -0.5*log(2*pi) + 0.5*(determinant(fullOpt$hessian)$mod - log(abs(tmpOpt$hessian))) - tmpOpt$obj + fullOpt$obj
      exp(ans)
    },
    error = function(e) {
      return(0.0)
    },
    warning = function(w) {
      return(NULL)
    },
    finally = {

    }
    )
    return(out)
  }

  fun.val <- sapply(X = par.val[,1], FUN =  marg2, ...)
  marg.sp <- splinefun(x = par.val[,1], y = fun.val)
  nc.marg <- integrate(marg.sp, lower = lo[1], upper = up[1])$value

  log.ncost <- log(nc.marg) + log.ncost
  log.den <- log(marg2(fullOpt$par[1], ...)) + log.den

  # the conditional
  cond2 = function(x, ...) {
    out <- tryCatch({
      tt = -ff(c(fullOpt$par[1], x), ...) + fullOpt$objective
      exp(tt)
    },
    error = function(e) {
      return(0.0)
    },
    warning = function(w) {
      return(NULL)
    },
    finally = {
    }
    )
    return(out)
  }

  lcfun.val <- sapply(X = par.val[,2], FUN = cond2, ...)
  lcond.sp <- splinefun(x = par.val[,2], y = lcfun.val)
  nc.cond <- integrate(lcond.sp, lower = lo[2], upper = up[2])$value

  log.ncost <- log(nc.cond) + log.ncost
  log.den <- log(cond2(fullOpt$par[2], ...)) + log.den

  return(-fullOpt$obj - log.den + log.ncost)
}
