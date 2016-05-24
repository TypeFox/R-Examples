#' @importFrom stats nobs 
#' @importFrom Rcpp evalCpp 
#' @useDynLib CompGLM

#' @title Conway-Maxwell Poisson GLM Fitting Function
#' @author Jeffrey Pollock <jeffpollock9@@gmail.com>
#' 
#' @description A function in similar format to \code{glm} which provides a linear 
#' form regressing on the parameters lambda and mu.
#' 
#' @details A log link is used for regression of the model parameters \eqn{\lambda} and \eqn{\nu}, that is:
#' \deqn{\log(\lambda) = \beta X}{log(\lambda) = \beta X} 
#' \deqn{\log(\nu) = \zeta Y}{log(\nu) = \zeta Y} 
#' where: \eqn{\beta} is the vector of coefficients for the parameter \eqn{\lambda},
#' \eqn{\zeta} is the vector of coefficients for the parameter \eqn{\nu},
#' \eqn{X} is the model matrix for the parameter \eqn{\lambda}, and
#' \eqn{Y} is the model matrix for the parameter \eqn{\nu}.
#' 
#' The parameter vectors are calculated via maximum likelihood using the general optimisation function
#' \code{\link{optim}}. A Poisson model will be fit using \code{\link{glm.fit}} and (unless starting values
#' are supplied) the coefficients will be used as starting values for the parameter vector \eqn{\beta}.
#' 
#' Several S3 functions have been implemented for model analysis \code{\link{print}}, \code{\link{coef}}, 
#' \code{\link{extractAIC}}, \code{\link{logLik}}, \code{\link{predict}}, and \code{\link{summary}},
#' 
#' @param lamFormula an object of class \code{\link{formula}} which determines the form of
#' regression for the model parameter \eqn{\lambda}. An offset can also be added in the formula.
#' @param nuFormula an object of class \code{\link{formula}} which determines the form of
#' regression for the model parameter \eqn{\nu}. The default value is \code{NULL} meaning the formula is
#' intercept only. An offset can also be added in the formula.
#' @param data an optional \code{\link{data.frame}} containing the variables in the model. If not found in data, 
#' the variables are taken from \code{environment(lamFormula)}.
#' @param lamStart optional vector of starting values for the coefficients of the \eqn{\lambda} regression.
#' @param nuStart optional vector of starting values for the coefficients of the \eqn{\nu} regression.
#' @param sumTo an integer for the summation term in the density (default 100).
#' @param method optimisation method passed to \code{\link{optim}} (default "BFGS").
#' @param ... further arguments to be passed to \code{\link{optim}}.
#' 
#' @return An object of class 'Comp' which is a list with all the components needed for the relevant S3 class
#' methods.
#' 
#' @references A Flexible Regression Model for Count Data, by Sellers & Shmueli, 
#' \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=1127359}
#' 
#' @examples 
#' set.seed(1)
#' n <- 5000
#' x1 <- rnorm(n, -1.0, 0.5)
#' x2 <- rnorm(n, 1.0, 0.7)
#' x3 <- rnorm(n, 2.0, 0.4)
#' y <- rpois(n, exp(-0.5 + 0.3 * x1 + 0.8 * x2 + 0.2 * x3))
#' data <- data.frame(y, x1, x2, x3)
#' model <- glm.comp(y ~ ., data = data)
#' print(model)
#' summary(model)
#' coef(model)
#' head(predict(model))
#' AIC(model)
#' @export glm.comp
glm.comp <- function(lamFormula, nuFormula = NULL, data, lamStart = NULL, nuStart = NULL, 
		sumTo = 100L, method = "BFGS", ...) {
	
	call <- match.call()
	
	if (missing(data)) {
		data <- environment(lamFormula)
	}
	
	lamModelFrame <- model.frame(lamFormula, data)
	lamModelTerms <- attr(lamModelFrame, "terms")
	
	y <- model.response(lamModelFrame)
	nobs <- NROW(y)
	
	if(any(y + 10 > sumTo)) {
		stop("Response within 10 of 'sumTo' argument")
	}
	
	xLam <- model.matrix(lamModelTerms, lamModelFrame)
	nLamVariables <- NCOL(xLam)
	
	if (is.null(nuFormula)) {
		xNu <- matrix(1.0, nobs, 1L)
		colnames(xNu) <- "(Intercept)"
		nuModelTerms <- NULL
		nuOffset <- rep.int(0L, nobs)
	} else {
		nuModelFrame <- model.frame(nuFormula, data)
		nuModelTerms <- attr(nuModelFrame, "terms")
		
		xNu <- model.matrix(nuModelTerms, nuModelFrame)
		nuOffset <- as.vector(model.offset(nuModelFrame))
	}
	
	nNuVariables <- NCOL(xNu)
	
	lamOffset <- as.vector(model.offset(lamModelFrame))
	if (is.null(lamOffset)) {
		lamOffset <- rep.int(0L, nobs)
	}
	
	poissonModel <- glm.fit(xLam, y, rep.int(1L, nobs), family = poisson(), offset = lamOffset)
	coefs <- poissonModel$coefficients
	if (!poissonModel$converged) {
		stop("attempt to fit poisson glm with glm.fit failed")
	}
	if (any(is.na(coefs))) { 
		badVariables = paste(names(coefs)[is.na(coefs)], collapse = ", ")
		warning(sprintf("initial parameter estimates return NA from glm.fit, dropping variables: %s", 
						badVariables)) 
		xLam <- xLam[ , !is.na(coefs), drop = FALSE]
	}
	
	if (is.null(lamStart)) {
		lamStart <- coefs[!is.na(coefs)]
	}
	
	if (is.null(nuStart)) {
		nuStart <-  rep.int(0L, nNuVariables)
		names(nuStart) <- colnames(xNu)
	}
	
	start <- c(nuStart, lamStart)
	
	nuIndexes <- seq_len(nNuVariables)
	lamIndexes <- seq_len(nLamVariables) + max(nuIndexes)
	
	fmin <- function(par) {
		
		zeta <- par[nuIndexes]
		nu <- exp(xNu %*% zeta + nuOffset)
		
		beta <- par[lamIndexes]
		lam <- exp(xLam %*% beta + lamOffset)
		
		pr <- dcomp(y, lam, nu, sumTo)
		return(-sum(log(pr)))
	}
	
	gmin <- function(par) {
		
		out <- numeric(length(par))
		
		zeta <- par[nuIndexes]
		nu <- exp(xNu %*% zeta + nuOffset)
		
		beta <- par[lamIndexes]
		lam <- exp(xLam %*% beta + lamOffset)
		
		betaGrad <- t(xLam) %*% (y - Y(lam, nu, sumTo) / Z(lam, nu, sumTo))
		zetaGrad <- t(xNu) %*% ((-log(factorial(y)) + W(lam, nu, sumTo) / Z(lam, nu, sumTo)) * nu)
		
		out[lamIndexes] <- betaGrad
		out[nuIndexes] <- zetaGrad
		
		return(-out)
	}
	
	res <- optim(start, fmin, gmin, method = method, hessian = TRUE, ...)
	
	beta <- res$par[lamIndexes]
	zeta <- res$par[nuIndexes]
	niter <- c(f.evals = res$counts[1L], g.evals = res$counts[2L])
	df.residual <- nobs - (nNuVariables + nLamVariables)
	
	fit <- list(call = call, beta = beta, zeta = zeta, logLik = -res$value,
			terms = lamModelTerms, nuModelTerms = nuModelTerms, data = data, lamOffset = lamOffset, 
			nuOffset = nuOffset, convergence = res$convergence, niter = niter, hessian = res$hessian, 
			df.residual = df.residual, df.null = nobs - 1L, nobs = nobs
	)
	
	class(fit) <- "Comp"
	return(fit)
}
