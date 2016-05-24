#' Calculate confidence interval for log determinant of covariance matrices

#' @param cov a covariance matrix or a list of covariance matrices
#' @param n   sample size, or vector of sample sizes
#' @param conf  confidence level
#' @param method 
#' @param bias.adj logical, to include bias correction term

#' @reference Cai, T. T., Liang, T. & Zhou, H. H. (2016). Law of log determinant of sample
#'         covariance matrix and optimal estimation of differential enteropy
#'         for high-dimensional Gaussian distributions.   JMA, to appear.
#'         www.stat.yale.edu/~hz68/Covariance-Determinant.pdf

logdetCI <- function(cov, n, conf=0.95, method=2, bias.adj=TRUE) {
	p <- if (is.list(cov)) sapply(cov, nrow)
	     else nrow(cov)
	if (length(p) > 1)
		if (!all(p == p[1])) stop("all covariance matrices must be the same size")
		else p <- p[1]

	bias <- 0
	if (method==2) {  # corollary 2
		if (bias.adj) {
			for (k in 1:p) bias <- bias + log( 1 - k/n)
		}
		se <- sqrt( -2 * log( 1 - p/n) )
	}
	else { # corollary 1
		if (bias.adj) bias <- p * (p+1) / (2*n)
		se <- sqrt(2 * p / n)
	}
	
	logdet <- if (is.list(cov)) sapply(cov, function(x) log( det(x) ))
	     else log( det(cov) )
	z <- qnorm((1-conf)/2, lower.tail=FALSE)
	lower <- (logdet - bias) - z * se
	upper <- (logdet - bias) + z * se
	
	result <- data.frame(logdet=logdet, bias=bias, se=se, lower=lower, upper=upper)
	result
}

# tau function, Eq(5), used as bias correction
tau <- function (n, p) {
	res <- 0
	for (k in 1:p) res <- res + (digamma( (n - k + 1) / 2) - log(n/2))
	return(res)
}	
	
