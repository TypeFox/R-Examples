#' @title Surface tests
#'
#' @description
#' Calculates surface parameters a1 to a4, the stationary point, the principal axes, the eigenvectors and -values
#'
#' @details
#' No details so far.
#'
#' @export
#' @param x Either an RSA object (returned by the \code{RSA} function), or the coefficient for the X predictor
#' @param y Y coefficient
#' @param x2 X^2 coefficient
#' @param y2 Y^2 coefficient
#' @param xy XY interaction coefficient
#' @param b0 The intercept
#' @param SE In case that the coefficients are provided directly (as parameters x, y, x2, y2, xy), SE can provide the standard errors of these estimates. SE has to be a named vector with exactly five elements with the names of the coefficients, e.g.: \code{SE=c(x=.1, y=.2, x2=.1, y2=.5, xy=.3)}. SEs of all parameters have to provided, otherwise the function will print an error. In case standard errors \emph{and} the covariances (see below) \emph{and} df (see below) are provided, parametric confidence intervals for a1 to a4 are calculated.
#' @param COV Covariances between parameters. COV has to be a named vector with exactly four elements with the names of four specific covariances, e.g.: \code{COV=c(x_y=.1, x2_y2 = .2, x2_xy = .3, y2_xy = .4)}, where x_y is the covariance between x and y, and so on. All these covariances have to provided with exactly these names, otherwise the function will print an error.
#' @param df Degrees of freedom for the calculation of a1 to a4 confidence intervals. The df are the residual dfs of the model (df = n - estimated parameters). For the full polynomial model, this is n - 6 in a regular regression (the following parameters are estimated: Intercept, x, y, xy, x2, y2). \code{df} should be a single number.
#' @param model If x is an RSA object, this parameter specifies the model from which to extract the coefficients

#' @return
#' Returns surface parameters a1 to a4. If an RSA object or SE, COV and df are provided, also significance test and standard errors of a1 to a4 are reported.
#' The stationary point (X0, Y0, and Z0). 
#' First principal axis (PA) relative to the X-Y plane (p10 = intercept, p11 = slope), second PA (p20 = intercept, p21 = slope). 
#' M = eigenvectors, l = eigenvalues, L = lambda matrix
#' as1X to as4X: surface parameters of the PA, relative to X values
#' as1Y to as4Y: surface parameters of the PA, relative to Y values
#' PA1.curvX: quadratic component of the first PA, as seen from X axis
#' PA2.curvX: quadratic component of the second PA, as seen from X axis
#' PA1.curv: quadratic component of the first PA, after optimal coord transformation
#' PA2.curv: quadratic component of the second PA, after optimal coord transformation
#'
#' @references
#' Shanock, L. R., Baran, B. E., Gentry, W. A., Pattison, S. C., & Heggestad, E. D. (2010). Polynomial Regression with Response Surface Analysis: A Powerful Approach for Examining Moderation and Overcoming Limitations of Difference Scores. \emph{Journal of Business and Psychology, 25}, 543-554. doi:10.1007/s10869-010-9183-4
#'
#' @seealso \code{\link{RSA}}
#'
#' @examples
#' # get surface parameters from known parameters
#' # example from Shanock et al. (2010), p. 548, Table 2
#' RSA.ST(x=-.23, y=.77, x2=-.07, y2=-.10, xy=.27)
#'
#'
#' ## Compute standard errors and p values for surface parameters 
#' ## from external regression coefficients:
#' # standard errors for coefficients
#' SE <- c(x=.09, y=.09, x2=.07, y2=.07, xy=.11)
#' # covariances for specific coefficients:
#' COV <- c(x_y= -.000, x2_y2 = .001, x2_xy = -.003, y2_xy = -.004)
#' RSA.ST(x = .131, y = .382, x2 = .074, xy = .002, y2 = .039, SE=SE, COV=COV, df=181)
#'
#'
#' # Get surface parameters from a computed RSA object
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 2
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	SD <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- SD + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df, models="full")
#' RSA.ST(r1)




RSA.ST <- function(x=0, y=0, x2=0, xy=0, y2=0, b0=0, SE=NULL, COV=NULL, df=NULL, model="full"){
	if (!is.null(attr(x, "class"))) {
		if (attr(x, "class") == "RSA") {
			if (model %in% c("absunc", "absdiff")) stop("No surface parameters are computed for absolute difference models!")
			if (inspect(x$models[[model]], "converged") == FALSE) {
				warning("Model has not converged!")
				return(NULL)
			}
			fit <- x
			C <- coef(fit$models[[model]])
			if (fit$models[[model]]@Options$estimator != "DWLS") {
				b0 <- as.numeric(ifelse(is.na(C[paste0(fit$DV, "~1")]), b0, C[paste0(fit$DV, "~1")]))
			} else {
				# the threshold is the negative of the intercept ...
				b0 <- -as.numeric(ifelse(is.na(C[paste0(fit$DV, "|t1")]), b0, C[paste0(fit$DV, "|t1")]))
			}
			x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
			y <- as.numeric(ifelse(is.na(C["b2"]), 0, C["b2"]))
			x2 <- as.numeric(ifelse(is.na(C["b3"]), 0, C["b3"]))
			y2 <- as.numeric(ifelse(is.na(C["b5"]), 0, C["b5"]))
			xy <- as.numeric(ifelse(is.na(C["b4"]), 0, C["b4"]))
			if (fit$models[[model]]@Options$estimator != "DWLS") {
				SE <- inspect(fit$models[[model]], "se")$beta[fit$DV, -1]
			} else {
				SE <- inspect(fit$models[[model]], "se")$gamma
			}
			names(SE) <- c("x", "y", "x2", "xy", "y2")
			COV0 <- vcov(fit$models[[model]])
			COV <- c(x_y = COV0["b1", "b2"], x2_y2=COV0["b3", "b5"], x2_xy = COV0["b3", "b4"], y2_xy = COV0["b5", "b4"])
			df <- lavaan::nobs(fit$models[[model]]) - sum(inspect(fit$models[[model]], "free")$beta != 0) - 1
		} else {
			stop(paste("Unknown object of class", attr(x, "class")))
		}
	} else {
		fit <- NULL
	}
	
	## compute stationary point and principal axes
	# Calculate slopes of principal axes (you can do that even when the stationary point is undefined)
	p11 <- (y2 - x2 + sqrt(((x2 - y2)^2) + xy^2)) / xy
	p21 <- (y2 - x2 - sqrt(((x2 - y2)^2) + xy^2)) / xy
	PA1.curvX <- x2 - xy*p11 + y2*p11^2
	PA2.curvX <- x2 - xy*p21 + y2*p21^2
	PA1.curv  <- PA1.curvX / p11^2
	PA2.curv  <- PA2.curvX / p21^2

	# calculate eigenvectors
	b <- c(x, y)	# linear terms
	B <- matrix(c(x2, xy, xy, y2), ncol=2, nrow=2) * matrix(c(1, .5, .5, 1), ncol=2)	# matrix of higher order terms

	M <- eigen(B)$vectors	# normalized eigenvectors
	l <- eigen(B)$values	# eigenvalues
	L <- t(M)%*%B%*%M		# Lambda matrix

	# Calculate stationary point, matrix-style
	# The next calculations only work when the stationary point is defined
	SP0 <- NULL
	try({SP0 <- -.5*solve(B)%*%b}, silent=TRUE)	
	if (!is.null(SP0)) {
		
		# Stationary point
		X0 <- SP0[1]
		Y0 <- SP0[2]		
		Z0 <- b0 + sum(c(x*X0, y*Y0, x2*X0^2, xy*X0*Y0, y2*Y0^2))
	
		# Intercepts of principal axes
		p10 <- as.numeric(Y0 - p11*X0)
		p20 <- as.numeric(Y0 - p21*X0)
		
		C1 <- -p10/(p11+1)
		C2 <- -p20/(p21+1)
		
		# a*: linear and quadratic slope of the principal axes
		as1X <- x + p11*y + xy*p10 + 2*y2*p10*p11
		as2X <- x2 + xy*p11 + (p11^2)*y2
		as1Y <- x/p11 + y - (2*x2*p10)/p11^2 - (xy*p10)/p11
		as2Y <- x2/p11^2 + xy/p11 + y2
	
		# second principal axis
		as3X <- x + p21*y + xy*p20 + 2*y2*p20*p21
		as4X <- x2 + xy*p21 + (p21^2)*y2
		as3Y <- x/p21 + y - (2*x2*p20)/p21^2 - (xy*p20)/p21
		as4Y <- x2/p21^2 + xy/p21 + y2
	} else {
		X0 <- Y0 <- Z0 <- p10 <- p20 <- C1 <- C2 <- as1X <- as2X <- as3X <- as4X <- as1Y <- as2Y <- as3Y <- as4Y <- NA
	}
	

		
	a1 <- as.numeric(x+y)
	a2 <- as.numeric(x2+y2+xy)
	a3 <- as.numeric(x-y)
	a4 <- as.numeric(x2-xy+y2)

	
	## SEs of Surface parameters: if an RSA object is provided, just retrieve them from that
	if (!is.null(fit)) {
		p.all <- as.data.frame(parameterEstimates(fit$model[[model]]))
		SP <- round(p.all[grepl("^a\\d$", p.all$label), c("est", "se", "z", "pvalue")], 10)
		colnames(SP) <- c("estimate", "SE", "z.value", "p.value")
	} else if (!is.null(SE) & length(SE) >= 5 & !is.null(COV) & length(COV) >= 4 & !is.null(df)) {
		## if SEs are provided: calculate confidence intervals and p values for a1 to a4
		se.a1 <- (sqrt(SE["x"]^2 + SE["y"]^2 + 2*COV["x_y"]))
		t.a1 <- a1 / se.a1
		p.a1 <- 2 * pt(-abs(t.a1), df)
		#se.a1 <- 
    
		se.a2 <- (sqrt(SE["x2"]^2 + SE["y2"]^2 + SE["xy"]^2 + 2*COV["x2_y2"] + 2*COV["x2_xy"] + 2*COV["y2_xy"]))
		t.a2 <- a2 / se.a2
		p.a2 <- 2 * pt(-abs(t.a2), df)

		se.a3 <- (sqrt(SE["x"]^2 + SE["y"]^2 - 2*COV["x_y"]))
		t.a3 <- a3 / se.a3
		p.a3 <- 2 * pt(-abs(t.a3), df)

		se.a4 <- (sqrt(SE["x2"]^2 + SE["y2"]^2 + SE["xy"]^2 - 2*COV["x2_xy"] + 2*COV["y2_xy"] - 2*COV["x2_y2"]))
		t.a4 <- a4 / se.a4
		p.a4 <- 2 * pt(-abs(t.a4), df)
	
		SP <- data.frame(estimate=c(a1, a2, a3, a4), SE=c(se.a1, se.a2, se.a3, se.a4), t.value=c(t.a1, t.a2, t.a3, t.a4), p.value=c(p.a1, p.a2, p.a3, p.a4))
		} else {
			SP <- data.frame(estimate=c(a1, a2, a3, a4), t.value=rep(NA, 4), p.value=rep(NA, 4))
		}
	
	rownames(SP) <- paste0("a", 1:4)
	
	PA <- data.frame(estimate=c(as1X, as2X, as3X, as4X, as1Y, as2Y, as3Y, as4Y), SE=NA, t.value=NA, p.value=NA)
	rownames(PA) <- c(paste0("as", 1:4, "X"), paste0("as", 1:4, "Y"))
		
	res <- list(X0=X0, Y0=Y0, Z0=Z0, p11=p11, p10=p10, p21=p21, p20=p20, C1=C1, C2=C2, SP=SP, PA=PA, M=M, l=l, L=L, PA1.curvX=PA1.curvX, PA2.curvX=PA2.curvX, PA1.curv=PA1.curv, PA2.curv=PA2.curv)
	return(res)
}
