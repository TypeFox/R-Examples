
#' Interpret bivariate threshold exceedance models
#'
#' This is an adaptation of the \pkg{evir} package \code{\link[evir]{interpret.gpdbiv}} function.
#' \code{interpret.fbvpot} was adapted to deal with the output of a call to
#' \code{\link[evd]{fbvpot}} from the \pkg{evd} and to handle families other than the logistic distribution.
#' The likelihood derivation comes from expression 2.10 in Smith et al. (1997).
#' @importFrom stats pbeta
#' @importFrom evd pbvevd
#' @seealso \code{\link[evir]{interpret.gpdbiv}}
#' @author Leo Belzile, adapting original S code by Alexander McNeil
#' @export
#' @references Smith, Tawn and Coles (1997), Markov chain models for threshold exceedances. \emph{Biometrika},
#' \strong{84}(2), 249--268.
#' @param fitted the output of \code{\link[evd]{fbvpot}} or a list. See Details.
#' @param q a vector of quantiles to consider, on the data scale. Must be greater than the thresholds.
#' @param silent boolean; whether to print the interpretation of the result. Default to \code{FALSE}.
#' @return an invisible numeric vector containing marginal, joint and conditional exceedance probabilities.
#' @details The list \code{fitted} must contain
#' \itemize{
#' \item \code{model} a string; see \code{\link[evd]{bvevd}} for options
#' \item \code{param} a named vector containing the parameters of the \code{model}, as well as parameters
#' \code{scale1}, \code{shape1},\code{scale2} and \code{shape2}, corresponding to marginal GPD parameters.
#' \item \code{threshold} a vector of length 2 containing the two thresholds.
#' \item \code{pat} the proportion of observations above the corresponding \code{threshold}
#' }
ibvpot <- function (fitted, q, silent=FALSE) {
	#If input is an object resulting from a call to "fpot"
	if(all(c("bvpot","evd") %in% class(fitted)) && length(class(fitted)==2)){
		fitted$pat <- fitted$nat[1:2]/nrow(fitted$data)
	}
	if(any(is.null(fitted$model),is.null(fitted$pat), is.null(fitted$threshold),
				is.null(fitted$param))){
	stop("Invalid argument for `fitted'. Please provide a list with components
							`model',`threshold', `pat',`param' or the output of `evd::fvbpot'")
	}
	if(! all(c("shape1","scale1","shape2","scale2") %in% names(fitted$param))){
		stop("Invalid arguments for `param'. Missing marginal parameters.")
	}
	if(length(fitted$pat)!=2 || length(fitted$threshold)!=2 || !is.numeric(fitted$threshold)
		|| any(fitted$pat>1) || any(fitted$pat<0)){
		stop("Invalid input for `pat' or `threshold' vector.")
	}
	if(! fitted$model %in% c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix")){
		stop("Invalid model. See documentation in the evd package")
	}
	if(fitted$model %in% c("log","hr","neglog")){
		modelclass <- "A"
	} else if(fitted$model %in% c("alog","aneglog")){
		modelclass <- "B"
	} else if(fitted$model %in% c("bilog", "negbilog","ct","amix")){
	 modelclass <- "C"
	}
	bivparnames <- switch(modelclass,
		A="dep",
		B=c("dep","asy1","asy2"),
    C=c("alpha","beta")
		)
	if(! all(bivparnames  %in% names(fitted$param))){
		stop("Invalid arguments for `param'. Missing bivariate parameters.")
	}
	#Exponent measure for the postulated model
  if(fitted$model=="ct" && !is.na(fitted$param['rho']) ){
    Vfuncf <- function(q, fitted){
      rho <- fitted$param['rho']; a1 <- fitted$param['alpha']; a2 <- fitted$param['beta']
      Ups <- exp((lbeta(a1,a2+rho)+log(q[2]))/rho)/
        (exp((lbeta(a1,a2+rho)+log(q[2]))/rho)+exp((lbeta(a1+rho,a2)+log(q[1]))/rho))
      exp(pbeta(Ups, a2,a1+rho, log.p=T)-log(q[1])) + exp(pbeta(1-Ups, a1,a2+rho, log.p=T)-log(q[2]))
    }
  } else{
    Vfuncf <- function(q, fitted){
      f <- switch(fitted$model,
                  log=pbvevd(q=q,model="log", dep=fitted$param['dep'], mar1=c(1,1,1)),
                  alog=pbvevd(q=q,model="alog", dep=fitted$param['dep'],asy=fitted$param[c('asy1','asy2')], mar1=c(1,1,1)),
                  hr=pbvevd(q=q,model="hr", dep=fitted$param['dep'], mar1=c(1,1,1)),
                  neglog=pbvevd(q=q,model="neglog", dep=fitted$param['dep'], mar1=c(1,1,1)),
                  aneglog=pbvevd(q=q,model="aneglog", dep=fitted$param['dep'],asy=fitted$param[c('asy1','asy2')], mar1=c(1,1,1)),
                  bilog=pbvevd(q=q,model="bilog", alpha=fitted$param['alpha'], beta=fitted$param['beta'], mar1=c(1,1,1)),
                  negbilog=pbvevd(q=q,model="negbilog", alpha=fitted$param['alpha'], beta=fitted$param['beta'], mar1=c(1,1,1)),
                  ct=pbvevd(q=q,model="ct", alpha=fitted$param['alpha'], beta=fitted$param['beta'], mar1=c(1,1,1)),
                  amix=pbvevd(q=q,model="amix", alpha=fitted$param['alpha'], beta=fitted$param['beta'], mar1=c(1,1,1))
      )
      return(-log(f))
    }
  }
	#Approximate joint CDF of bivariate EV model: F(x,y)
  bivcdf <- function(x, y,  u1, lambda1,  sigma1, xi1,
                      u2, lambda2, sigma2, xi2,  vfunc){
	    Zfunc <- function(y, u, lambda, xi, sigma){
	    	(lambda^-1) * exp(log(1 + (xi * pmax((y - u), 0))/sigma)/xi)
	    }
    1 - vfunc(c(Zfunc(x, u1, lambda1, xi1, sigma1),
    						Zfunc(y, u2, lambda2, xi2, sigma2)), fitted)
  }
  marg <- function(x, u1, lambda1, sigma1,xi1) {
    1 - lambda1 * exp(-log(1 + (xi1 * (x - u1))/sigma1)/xi1)
  }
  bivsurv <- function(x, y,  u1, lambda1, xi1, sigma1,
                       u2, lambda2, xi2, sigma2, marg, newfunc, vfunc) {
    1 - marg(x, u1, lambda1, xi1, sigma1) - marg(y, u2, lambda2, xi2, sigma2) +
  		newfunc(x, y, u1, lambda1, xi1, sigma1, u2, lambda2, xi2, sigma2, vfunc)
  }
  if (fitted$threshold[1] > q[1])
    stop("Point below x threshold")
  if (fitted$threshold[2] > q[2])
    stop("Point below y threshold")
  p1 <- 1 - marg(q[1], fitted$threshold[1], fitted$pat[1], fitted$param['scale1'],
                 fitted$param['shape1']) #u, lambda, scale, shape
  p2 <- 1 - marg(q[2], fitted$threshold[2], fitted$pat[2], fitted$param['scale2'],
                 fitted$param['shape2'])
  p12 <- bivsurv(q[1], q[2],  fitted$threshold[1], fitted$pat[1], fitted$param['scale1'],
                  fitted$param['shape1'], fitted$threshold[2], fitted$pat[2], fitted$param['scale2'],
                  fitted$param['shape2'], marg, bivcdf, Vfuncf)
  if(!silent){
  cat("Bivariate POT model:", fitted$model, "\n")
  cat("Thresholds:", fitted$threshold[1], fitted$threshold[2], "\n")
  cat("Extreme levels of interest (x,y):", q[1], q[2], "\n")
  cat("P(X > x) =", format(p1), "\n")
  cat("P(Y > y) =", format(p2), "\n")
  cat("P(X > x, Y > y) =", format(p12), "\n")
  cat("P(X > x) * P(Y > y) =", format(p1 * p2), "\n")
  cat("P(Y > y | X > x) =", format(p12/p1), "\n")
  cat("P(X > x | Y > y) =", format(p12/p2), "\n")
  }
  output <- as.numeric(c(p1, p2, p12, p1*p2, p12/p1, p12/p2))
  names(output) <- c("p(1)","p(2)","p(1,2)","p(1)p(2)","p(2|1)", "p(1|2)")
  invisible(output)

}
