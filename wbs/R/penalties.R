#' @rdname ssic.penalty
#' @title Strengthened Schwarz Information Criterion penalty term
#' @description The function evaluates the penalty term for the strengthened Schwarz Information Criterion proposed in P. Fryzlewicz (2014). This routine is typically not called directly by the user; its name can be passed as an argument to \code{\link{changepoints}}.
#' @param n the number of observations
#' @param cpt a vector with localisations of change-points
#' @param alpha a scalar greater than one
#' @param ssic.type a string ("log" or "power")
#' @export ssic.penalty
#' @return the penalty term \eqn{k(\log(n))^{alpha}}{k(log(n))^(alpha)} for \code{ssic.penalty="log"} or \eqn{k n^{alpha}}{k * n^(alpha)} for \code{ssic.penalty="power"}, where \eqn{k}{k} denotes the number of elements in \code{cpt}
#' @references P. Fryzlewicz (2014), Wild Binary Segmentation for multiple change-point detection. Annals of Statistics, to appear. (\url{http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf})
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- wbs(x)
#' w.cpt <- changepoints(w,penalty="ssic.penalty")
#' w.cpt$cpt.ic

ssic.penalty <- function(n,cpt,alpha=1.01,ssic.type=c("log","power")){
	alpha <- as.numeric(alpha)

	ssic.type <- match.arg(ssic.type)
	
	if (ssic.type == "log")  pen <- log(n)^alpha
	else if (ssic.type == "power") n^alpha
	
	k <- length(cpt)
	
	return(k*pen)
}


#' @rdname bic.penalty
#' @title Bayesian Information Criterion penalty term
#' @description  The function evaluates the penalty term for the standard Bayesian Information Criterion applied to the change-point detection problem. This routine is typically not called directly by the user; its name can be passed as an argument to \code{\link{changepoints}}.  
#' @return the penalty term \eqn{k\log(n)}{k * log(n)} where  \eqn{k}{k} denotes the number of elements in \code{cpt}
#' @param n the number of observations
#' @param cpt a vector with localisations of change-points
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- wbs(x)
#' w.cpt <- changepoints(w,penalty="bic.penalty")
#' w.cpt$cpt.ic
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- wbs(x)
#' w.cpt <- changepoints(w,penalty="bic.penalty")
#' w.cpt$cpt.ic

bic.penalty <- function(n,cpt){
	k <- length(cpt)
	return(k*log(n))
}

#' @rdname mbic.penalty
#' @title Modified Bayes Information Criterion penalty term
#' @description  The function evaluates the penalty term for the  Modified Bayes Information Criterion proposed in N. Zhang and D. Siegmund (2007). This routine is typically not called directly by the user; its name can be passed as an argument to \code{\link{changepoints}}.
#' @references N. Zhang and D. Siegmund (2007), A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data, Biometrics.
#' @param n the number of observations
#' @param cpt a vector with localisations of change-points
#' @return the penalty term \deqn{\frac{3}{2}k\log(n)+\frac{1}{2}\sum_{i=1}^{k+1}\log\frac{l_{i}}{n},}{3/2 * k * log(n)+1/2 * sum_i^k+1 log(l_i)/n,} where  \eqn{k}{k} denotes the number of elements in \code{cpt} and \eqn{l_{i}}{l_i} are the lengths of the intervals between changepoints in \code{cpt}
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- wbs(x)
#' w.cpt <- changepoints(w,penalty="mbic.penalty")
#' w.cpt$cpt.ic

mbic.penalty <- function(n,cpt){
	k <- length(cpt)
	if(k>0)	3/2 * k * log(n) + 1/2 * sum(log(diff(c(0,sort(cpt),n))/n))
	else 1/2 * sum(log(diff(c(0,n))/n))
}

