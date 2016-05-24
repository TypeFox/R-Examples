# function for doing multivariate non-linear tests on parameters
# tests the probability that fun(beta) > 0, where beta are the estimated parameters
# from felm()






#' Compute expectation of a function of the coefficients.
#' 
#' Use integration of the joint distribution of the coefficients to compute the
#' expectation of some function of the coefficients.  Can be used for
#' non-linear inference tests.
#' 
#' The function \code{nlexpect} integrates the function \code{fun(x)} over the
#' multivariate normal distribution specified by the point estimates and the
#' covariance matrix \code{vcov(est)}.  This is the expectation of
#' \code{fun(beta)} if we were to bootstrap the data (e.g. by drawing the
#' residuals anew) and do repeated estimations.
#' 
#' The list of coefficients used by \code{fun} must be specified in
#' \code{coefs}.
#' 
#' If the function is simple, it can be specified as a quoted expression like
#' \code{quote(a*b+log(abs(d)))}. In this case, if \code{coefs} is not
#' specified, it will be set to the list of all the variables occuring in the
#' expression which are also names of coefficients.
#' 
#' \code{fun} may return a vector of values, in which case a vector of
#' expectations is computed, like \code{quote(c(a*b, a^3-b))}. However, if the
#' expressions contain different variables, like \code{quote(c(a*b, d*e))}, a
#' quite compute intensive 4-dimensional integral will be computed, compared to
#' two cheap 2-dimensional integrals if you do them separately.
#' 
#' You may of course also integrate inequalites like \code{quote(abs(x1-0.2) >
#' 0.2)} to simulate the probability from t-tests or Wald-tests. See the
#' examples.
#' 
#' The function you provide will get an argument \code{...} if it does not have
#' one already.  It will also be passed an argument \code{.z} which contains
#' the actual coefficients in normalized coordinates, i.e. if \code{ch} is the
#' Cholesky decomposition of the covariance matrix, and \code{pt} are the point
#' estimates, the coefficients will be \code{pt + ch \%*\% .s}
#' 
#' The \code{tol} argument specifies both the relative tolerance and the
#' absolute tolerance. If these should not be the same, specify \code{tol} as a
#' vector of length 2. The first value is the relative tolerance, the second is
#' the absolute tolerance. Termination occurs when at least one of the
#' tolerances is met.
#' 
#' The \code{...} can be used for passing other arguments to the integration
#' routine \code{\link[R2Cuba]{cuhre}}.
#' 
#' @param est object of class \code{"felm"}, a result of a call to
#' \code{\link{felm}}.
#' @param fun function of coefficients to be integrated. Can also be a
#' \code{quote}d expression.
#' @param coefs character. Names of coefficients to test. Only needed if
#' \code{fun} is a function.
#' @param ... other arguments passed to fun or the integration routine.
#' @param tol numeric. Tolerance for the computed probability.
#' @param lhs character. Name of the left hand side, if \code{est} has more
#' than one.
#' @param cv Covariance matrix to use in place of \code{vcov(est)}
#' @param istats logical. Should convergence information from the integration
#' routine be included as attributes?
#' @param flags list. Additional flags for the underlying integration routine.
#' @param max.eval integer. Maximum number of integral evaluations.
#' @return The function \code{nlexpect} computes and returns the expectation of
#' the function \code{fun(beta)}, with \code{beta} a vector of coefficients.
#' I.e., if the coefficients \code{beta} are bootstrapped a large number of
#' times, \code{nlexpect(est, fun)} should be equal to \code{mean(fun(beta))}.
#' Additional diagnostic output from \code{\link[R2Cuba]{cuhre}} is returned in
#' attributes if \code{istats=TRUE}.
#' @note An alternative to this method is to use the \code{bootexpr} argument
#' with \code{\link{felm}}, to do a Monte Carlo integration.
#' 
#' This function requires the package \pkg{R2Cuba}, and will fail if it is not
#' available.
#' @seealso \code{\link{waldtest}}
#' @examples
#' 
#' N <- 100
#' x1 <- rnorm(N)
#' # make some correlation
#' x2 <- 0.1*rnorm(N) + 0.1*x1
#' y <- 0.1*x1 + x2 + rnorm(N)
#' summary(est <- felm(y ~ x1 + x2))
#' pt1 <- coef(est)['x1']
#' pt2 <- coef(est)['x2']
#' # expected values of coefficients, should match the summary
#' # and variance, i.e. square of standard errors in the summary
#' nlexpect(est, quote(c(x1=x1,x2=x2,var=c((x1-pt1)^2,(x2-pt2)^2))))
#' \donttest{
#' # the covariance matrix:
#' nlexpect(est, tcrossprod(as.matrix(c(x1-pt1,x2-pt2))))
#' }
#' #Wald test of single variable
#' waldtest(est, ~x1)['p.F']
#' # the same with nlexpect, i.e. probability for observing abs(x1)>abs(pt1) conditional
#' # on E(x1) = 0.
#' nlexpect(est, (x1-pt1)^2 > pt1^2)
#' # which of course is the same as
#' 2*nlexpect(est, x1 < 0)
#' 
#' \donttest{
#' # then a joint test. Here we find the probability that
#' # we are further from origo than the point estimates, measured
#' # in a variance normalized coordinate system.
#' waldtest(est, ~ x1 | x2)['p.F']
#' #inverse cholesky
#' ich <- solve(t(chol(vcov(est)[c('x1','x2'), c('x1','x2')])))
#' # convert to normalized coordinates
#' nz <- sum((ich %*% c(pt1,pt2))^2)
#' # find probability that we're further away
#' nlexpect(est, sum((ich %*% c(x1-pt1,x2-pt2))^2) > nz)
#' # or use the .z argument provided automatically
#' nlexpect(est, sum(.z^2) > nz, coefs=c('x1','x2'))
#' 
#' # Non-linear test:
#' f <- function(x) c(poly=x[['x1']]*(6*x[['x1']]-x[['x2']]^2))
#' waldtest(est, f)['p.F']
#' # In general, for a function f, the non-linear Wald test is something like
#' # the following:
#' # expected value of function
#' Ef <- nlexpect(est, f, coefs=c('x1','x2'))
#' # point value of function
#' Pf <- f(c(pt1,pt2))
#' # similar to a Wald test:
#' nlexpect(est, function(x) (f(x)-Ef)^2 > Pf^2, c('x1','x2'))
#' # one-sided
#' nlexpect(est, function(x) f(x)-Ef > abs(Pf), c('x1','x2'))
#' # other sided
#' nlexpect(est, function(x) f(x)-Ef < -abs(Pf), c('x1','x2'))
#' }
#' 
#' @export nlexpect
nlexpect <- function(est, fun, coefs, ..., tol=getOption('lfe.etol'), lhs=NULL,
                     cv,istats=FALSE,flags=list(verbose=0), max.eval=100000L) {
#  if(!requireNamespace('cubature', quietly=TRUE)) {warning('Package "cubature" not found.'); return(NULL);}
  if(!requireNamespace('R2Cuba', quietly=TRUE)) {warning('Package "R2Cuba" not found.'); return(NULL);}
  if(isTRUE(est$nostats) && missing(cv))
      stop('This test requires that felm() is run without nostats=TRUE; or specify a cv argument')
  # Find the covariance matrix
  if(missing(cv)) cv <- vcov(est, lhs=lhs)

  # Some kludge to be able to use non-quoted expression for fun
  afun <- substitute(fun)
  if(is.call(afun) || (is.name(afun) && (as.character(afun) %in% colnames(cv)))) {
    lfun <- as.list(afun)
    if(identical(lfun[[1]], quote(expression)))
        fun <- as.call(lfun[[2]])
    else if(!identical(lfun[[1]],quote(quote)) && !identical(lfun[[1]],quote(`function`)))
        fun <- afun
  }

  if(is.call(fun) || is.name(fun)) {
    # it's an expression. Figure out the coefficients used
    if(missing(coefs)) coefs <- intersect(all.vars(fun),colnames(cv))
    # make it a function
    fun <- local(function(x, ...) eval(fun,c(as.list(x),list(...))), list(fun=fun))
  } else if(is.function(fun)) {
    #add a ... formal if it doesn't exist
    fa <- formals(fun)
    if(!('.z' %in% names(fa)))
        formals(fun) <- c(fa,alist(.z=))
    if(!('...' %in% names(fa)))
        formals(fun) <- c(fa,alist(...=))
  }

  if(missing(coefs) || length(coefs %in% colnames(cv)) == 0)
      stop('No coefficients specified')
  # Find the coefficients
  cf <- drop(coef(est, lhs=lhs))[coefs]
  # and the Cholesky
  ch <- t(chol(cv[coefs,coefs,drop=FALSE]))

  # Now, we need to integrate fun(x) > 0 over the joint distribution of the parameters
  # We do this as follows. We integrate over a standard hypercube (-1,1) x (-1,1) x ...
  # adaptIntegrate can't take infinite limits.
  # We first transform these to (-Inf, Inf) with
  # z = x/(1-x^2)
  # the Jacobian determinant becomes the product of 
  # (1+z^2)/((1-z^2)^2)
  # We transform the integration variables with the covariance matrix to feed fun(),
  # then integrate fun(x) > 0 with the multivariate normal distribution.
  # we use the package cubature for the integration.
  
  integrand <- function(x, ...) {
    x2 <- x^2L
    jac <- prod((1+x2)/((1-x2)^2L))
    z <- x/(1-x2)
    # z is the standard normal (t really) multivariate
#    dens <- prod(dnorm(z))
    dens <- prod(dt(z,est$df))
    beta <- drop(cf + ch %*% z)
    names(beta) <- coefs
    ret <- fun(beta, .z=z, ...)*jac*dens
    if(anyNA(ret)) stop('Function value is NA for beta: ',format(beta),' ',x)
    ret
  }

  K <- length(cf)
  sv <- fun(cf, .z=rep(0,K), ...)
  fdim <- length(sv)
  if(length(tol) == 2) {
    reltol <- tol[1]
    abstol <- tol[2]
  } else {
    reltol <- abstol <- tol
  }
  ## ret <- cubature::adaptIntegrate(integrand,rep(-1,K),rep(1,K),...,tol=reltol,absError=abstol,fDim=fdim)
  ## names(ret$integral) <- names(sv)
  ## dim(ret$integral) <- dim(sv)
  ## if(is.array(sv)) dimnames(ret$integral) <- dimnames(sv)
  ## if(!istats) return(ret$integral)
  ## names(ret)[match('integral',names(ret))] <- names(as.list(args(structure)))[1]  
  ret <- R2Cuba::cuhre(K,fdim,integrand,lower=rep(-1,K),...,
                       upper=rep(1,K),rel.tol=reltol,abs.tol=abstol,
                       flags=flags, max.eval=max.eval)
  names(ret$value) <- names(sv)
  if(is.array(sv)) {
    dim(ret$value) <- dim(sv)
    dimnames(ret$value) <- dimnames(sv)
  }
  if(ret$ifail != 0) warning('integration failed with: "',ret$message, '", use istats=TRUE to see details')
  if(!istats) return(ret$value)
  names(ret)[match('value',names(ret))] <- names(as.list(args(structure)))[1]
  do.call(structure,ret)
}
