calcBetaParam <- function(mode, strength) {
  a <- plogis(mode) * strength
  b <- strength - a
  c(a=a, b=b, c=log(beta(a+1,b+1)))
}

#' Univariate priors commonly used in IFA models
#'
#' The returned model evaluates to the fit of the priors in deviance
#' (-2 log likelihood) units. The analytic gradient and Hessian are
#' included for quick optimization using Newton-Raphson.
#'
#' Priors of type 'beta' and 'logit-norm' are commonly used for the
#' lower asymptote parameter of the 3PL model. Both of these priors
#' assume that the parameter is in logit units. The 'lnorm' prior
#' can be used for slope parameters.
#'
#' @param type one of c("lnorm","beta","logit-norm")
#' @param labels a vector of parameters to which to apply the prior density
#' @param mode the mode of the prior density
#' @param strength a prior-specific strength (optional)
#' @param name the name of the mxModel returned
#' @return
#' an mxModel that evaluates to the prior density in deviance units
#' @export
#' @import OpenMx
#' @importFrom stats plogis
#' @examples
#' model <- univariatePrior("logit-norm", "x1", -1)
#' model$priorParam$values[1,1] <- -.6
#' model <- mxRun(model)
#' model$output$fit
#' model$output$gradient
#' model$output$hessian
univariatePrior <- function(type, labels, mode, strength=NULL, name="univariatePrior") {
	priorParam <- mxMatrix(name='priorParam', nrow=1, ncol=length(labels),
			       free=TRUE, labels=labels)
	priorParam$values[,] <- NA
	if (type == "beta") {
		if (length(strength) == 0) {
			strength <- 5
		}
		betaParam <- mapply(calcBetaParam, mode, strength)
		betaA <- mxMatrix(name='betaA', values=betaParam['a',,drop=FALSE])
		betaB <- mxMatrix(name='betaB', values=betaParam['b',,drop=FALSE])
		betaC <- mxMatrix(name='betaC', values=betaParam['c',,drop=FALSE])
		betaFit <- mxAlgebra(2 * sum((betaA + betaB)*log(exp(priorParam)+1) -
					     betaA * priorParam + betaC), name='betaFit')
		betaGrad <- mxAlgebra(2*(betaB-(betaA+betaB)/(exp(priorParam) + 1)),
				      name='betaGrad', dimnames=list(c(),priorParam$labels))
		betaHess <- mxAlgebra(vec2diag(2*(exp(priorParam)*(betaA+betaB) / (exp(priorParam) + 1)^2)),
				      name='betaHess',
				      dimnames=list(priorParam$labels, priorParam$labels))
		mxModel(model=name, priorParam, betaA, betaB, betaC,
			betaFit, betaGrad, betaHess,
			mxFitFunctionAlgebra('betaFit', gradient='betaGrad', hessian='betaHess'))
	} else if (type == "logit-norm") {
		if (length(strength) == 0) {
			strength <- .5
		}
		gaussM <- mxMatrix(name='gaussM', nrow=1, ncol=length(mode), values=mode)
		gaussSD <- mxMatrix(name='gaussSD', nrow=1, ncol=length(mode), values=strength)
		gaussFit <- mxAlgebra(sum(log(2*pi) + 2*log(gaussSD) +
					  (priorParam-gaussM)^2/gaussSD^2), name='gaussFit')
		gaussGrad <- mxAlgebra(2*(priorParam - gaussM)/gaussSD^2, name='gaussGrad',
				       dimnames=list(c(),priorParam$labels))
		gaussHess <- mxAlgebra(vec2diag(2/gaussSD^2), name='gaussHess',
				       dimnames=list(priorParam$labels, priorParam$labels))
		mxModel(model=name, priorParam, gaussM, gaussSD,
			gaussFit, gaussGrad, gaussHess,
			mxFitFunctionAlgebra('gaussFit', gradient='gaussGrad', hessian='gaussHess'))
	} else if (type == "lnorm") {
		if (length(strength) == 0) {
			strength <- 1
		}
		meanlog <- mxMatrix(name='meanlog', nrow=1, ncol=length(mode), values=log(mode) + strength^2)
		sdlog <- mxMatrix(name='sdlog', nrow=1, ncol=length(mode), values=strength)
		fit <- mxAlgebra(sum(log(2*pi) + 2*log(sdlog) + 2*log(priorParam) +
		                       (log(priorParam) - meanlog)^2/(sdlog*sdlog)), name="fit")
		grad <- mxAlgebra(2*(log(priorParam) - meanlog)/(sdlog*sdlog * priorParam) +
		                    2/priorParam, name='grad',
		                  dimnames=list(c(),priorParam$labels))
		hess <- mxAlgebra(vec2diag(-2*(log(priorParam) - meanlog)/(sdlog*sdlog * priorParam*priorParam) +
		                             2/(sdlog*sdlog * priorParam*priorParam) - 2/(priorParam*priorParam)),
		                  name='hess',
		                  dimnames=list(priorParam$labels, priorParam$labels))
		mxModel(model=name, priorParam, meanlog, sdlog, fit, grad, hess,
			mxFitFunctionAlgebra('fit', gradient='grad', hessian='hess'))
	} else {
		stop(paste("Unknown prior type", type))
	}
}

#' Uniqueness prior to assist in item factor analysis
#' 
#' To prevent Heywood cases, Bock, Gibbons, & Muraki (1988)
#' suggested a beta prior on the uniqueness (Equations 43-46).
#' The analytic gradient and Hessian are
#' included for quick optimization using Newton-Raphson.
#' 
#' To reproduce these derivatives in maxima for the case
#' of 2 slopes (\code{c} and \code{d}), use the following code:
#' 
#' \code{f(c,d) := -p*log(1-(c^2 / (c^2+d^2+1) + (d^2 / (c^2+d^2+1))));}
#' 
#' \code{diff(f(c,d), d),radcan;}
#' 
#' \code{diff(diff(f(c,d), d),d),radcan;}
#' 
#' The general pattern is given in Bock, Gibbons, & Muraki.
#' 
#' @param model an mxModel
#' @param numFactors the number of factors.
#' All items are assumed to have the same number of factors.
#' @param strength the strength of the prior
#' @param name the name of the mxModel that is returned
#' @return an mxModel that evaluates to the prior density in deviance units
#' @export
#' @references
#' Bock, R. D., Gibbons, R., & Muraki, E. (1988). Full-information item factor analysis.
#' \emph{Applied Psychological Measurement, 12}(3), 261-280.
#' @examples
#' numItems <- 6
#' spec <- list()
#' spec[1:numItems] <- rpf.drm(factors=2)
#' names(spec) <- paste0("i", 1:numItems)
#' item <- mxMatrix(name="item", free=TRUE,
#'                  values=mxSimplify2Array(lapply(spec, rpf.rparam)))
#' item$labels[1:2,] <- paste0('p',1:(numItems * 2))
#' data <- rpf.sample(500, spec, item$values)
#' m1 <- mxModel(model="m1", item,
#'               mxData(observed=data, type="raw"),
#'               mxExpectationBA81(spec),
#'               mxFitFunctionML())
#' up <- uniquenessPrior(m1, 2)
#' container <- mxModel("container", m1, up,
#'                      mxFitFunctionMultigroup(c("m1", "uniquenessPrior")),
#'                      mxComputeSequence(list(
#'                        mxComputeOnce('fitfunction', c('fit','gradient')),
#'                        mxComputeReportDeriv())))
#' container <- mxRun(container)
#' container$output$fit
#' container$output$gradient
uniquenessPrior <- function(model, numFactors, strength=0.1,
                            name="uniquenessPrior") {
  ex <- model$expectation
  if (!is(ex, "MxExpectationBA81")) {
    stop(paste(model$name, "does not contain MxExpectationBA81"))
  }

  imat <- model[[ex$item]]
  if (is.null(imat)) {
    stop(paste("Cannot find item matrix in model", model$name))
  }

  numItems <- ncol(imat)
  imatName <- paste(model$name, imat$name, sep=".")
  str <- c("1")
  for (fx in 1:numFactors) {
    term <- paste(imatName, "[", fx, ", 1:",  numItems, "]")
    str <- c(str, " + ", term, "*" ,term)
  }
	den <- mxAlgebraFromString(paste(str, collapse=""), name='den')
	den2 <- mxAlgebra(den*den, name='den2')
	
	str <- c("-", strength, "*(0")
	for (ix in 1:numItems) {
	  term <- paste(imatName,"[1:",numFactors,",",ix,"]", sep="")
	  str <- c(str,"+log(1-sum(",term,"*",term,"/den[1,",ix,"]))")
	}
	str <- c(str, ")")
	fit <- mxAlgebraFromString(paste(str, collapse=""), name="fit")
	
	mask <- matrix(FALSE, nrow=nrow(imat), ncol=ncol(imat))
	mask[1:numFactors,] <- imat$free[1:numFactors,]
	fv <- unique(imat$labels[mask])
	if (any(is.na(fv))) stop("All free slope parameters must be labelled")
	if (!length(fv)) stop("No free parameters")
	
	fvMap <- mxMatrix(values=0, nrow=numFactors * numItems, ncol=length(fv), name="fvMap")

	mx1 <- 1
	str1 <- paste(strength," * 2 * cbind(", sep="")
	str2 <- paste("2 * ",strength,"*rbind(", sep="")
	for (ix1 in 1:numItems) {
	  for (fx1 in 1:numFactors) {
	    str2 <- paste(str2, "cbind(", sep="")
	    if (!imat$free[fx1,ix1]) {
	      str1 <- paste(str1, 0)
	    } else {
	      str1 <- paste(str1, imatName,"[",fx1,",",ix1,"]/den[1,",ix1,"]", sep="")
	      fvMap$values[mx1, match(imat$labels[fx1,ix1], fv)] <- 1
	    }
	    for (ix2 in 1:numItems) {
	      for (fx2 in 1:numFactors) {
	        if (!imat$free[fx2,ix2] || ix1 != ix2) {
	          str2 <- paste(str2, 0)
	        } else {
	          if (fx1==fx2) {
	            term <- paste(imatName,"[",fx1,",",ix1,"]", sep="")
	            str2 <- paste(str2, "(den[1,",ix1,"] - 2*",term,"*",term,")/den2[1,",ix1,"]",sep="")
	          } else {
	            term1 <- paste(imatName,"[",fx1,",",ix1,"]", sep="")
	            term2 <- paste(imatName,"[",fx2,",",ix2,"]", sep="")
	            str2 <- paste(str2, "- 2*",term1,"*",term2,"/den2[1,",ix1,"]",sep="")
	          }
	        }
	        if (ix2*fx2 < numItems * numFactors) {
	          str2 <- paste(str2, ", ", sep="")
	        }
	      }
	    }
	    str2 <- paste(str2, ")", sep="")
	    if (mx1 < numItems * numFactors) {
	      str1 <- paste(str1, ", ", sep="")
	      str2 <- paste(str2, ", ", sep="")
	    }
	    mx1 <- mx1+1
	  }
	}
	str1 <- paste(str1, ")", sep="")
	str2 <- paste(str2, ")", sep="")

	allGrad <- 	mxAlgebraFromString(str1, name='allGrad')
	grad <- mxAlgebra(allGrad %*% fvMap, name='grad', dimnames=list(c(), fv))
	allHess <- mxAlgebraFromString(str2, name='allHess')
	hess <- mxAlgebra(t(fvMap) %*% allHess %*% fvMap, name='hess', dimnames=list(fv, fv))
	
	mxModel(model=name, den, den2, fit, fvMap, allGrad, grad, allHess, hess,
	        mxFitFunctionAlgebra('fit', gradient='grad', hessian='hess'))
}
