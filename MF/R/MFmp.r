
#' @description Estimates mitigated fraction from matched pairs.
#' @details Estimates \emph{MF} from matched pairs by the difference of multinomial fractions \eqn{(\Sigma I(x<y) - \Sigma I(x>y)) / N}. The trinomial vector is \eqn{\{\Sigma I(x<y), \Sigma I(x=y), \Sigma I(x>y)\}}
#' @title Mitigated fraction from matched pairs
#' @usage MFmp(formula=NULL, data=NULL, compare = c("con", "vac"), x=NULL, alpha=0.05, df=NULL, tdist=T)
#' @param formula Formula of the form \code{y ~ x + cluster(w)}, where y is a continuous response, x is a factor with two levels of treatment, and w is a factor indicating the clusters.
#' @param data Data frame
#' @param compare Text vector stating the factor levels - \code{compare[1]} is the control or reference group to which \code{compare[2]} is compared
#' @param x Trinomial vector \eqn{\{\Sigma I(x<y), \Sigma I(x=y), \Sigma I(x>y)\}}
#' @param alpha Complement of the confidence level.
#' @param df Degrees of freedom. Default N-2
#' @param tdist Use quantiles of t or Gaussian distribution for confidence interval? Default t distribution.
#' @export
#' @return a \code{\link{mfmp-class}} data object
#' @seealso \code{\link{mfmp-class}}
#' @references Siev D. (2005). An estimator of intervention effect on disease severity. \emph{Journal of Modern Applied Statistical Methods.} \bold{4:500--508}
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @examples
#' MFmp(les ~ tx + cluster(cage), mlesions, compare = c('con', 'vac'))
#' MFmp(x = c(12, 12, 2))
MFmp <- function(formula = NULL, data = NULL, compare = c("con", "vac"), x = NULL, 
	alpha = 0.05, df = NULL, tdist = T){
	# asymptotic CI for matched pairs
	# x is a trinomial frequency vector
	# c(x>y,x=y,x<y))
	# difference of multinomial fractions
	# I(x<y) - I(x>y)

	if(!is.null(formula) & !is.null(data)){
		byCluster <- MFClus(formula = formula, data = data, compare = compare)$byCluster[, 'mf']
		x <- rev(c(table(byCluster)))
	} else if(is.null(x)) {
		stop('Need to supply either formula and data or x = vector')
	}
	
	N <- sum(x)
	p <- x/N
	V <- (diag(p) - t(t(p)) %*% t(p)) / N
	V 
	A <- grad <- c(1, 0, -1)
	B <- t(A) %*% p
	VB <- t(A) %*% V %*% A

	gradl <- c(1 / (p[1] - p[3]), 0, 1 / (p[1] - p[3]))
	logB <- log(B)
	VlogB <- t(gradl) %*% V %*% gradl

	if(tdist & is.null(df)){
		df <- N - 2
	}
	
	if(!is.null(df)){
		q <- qt(c(0.5, alpha/2, 1 - alpha/2), df)
		what <- paste(100 * (1 - alpha), "% t intervals on ", df, " df\n", sep = "")
	} else {
		q <- qnorm(c(0.5, alpha/2, 1 - alpha/2))
		what <- paste(100 * (1 - alpha), "% gaussian interval\n", sep = "")
	}

	ci <- B + q * sqrt(VB) 
	names(ci) <- c("point", "lower", "upper")

	# out <- list(ci = ci, x = x, what = what, alpha = alpha, tdist = tdist, df = df)
	# class(out) <- 'mfmp'
	# return(out)
	return(mfmp$new(ci = ci, x = x, what = what, alpha = alpha, tdist = tdist, df = df))
}
	
