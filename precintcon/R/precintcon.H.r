#' @noRd
#' @name precintcon.H
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.H 
#' @title The cumulative probability.
#' @description The cumulative probability defined by
#' \deqn{pzero + (1 - pzero) \times pgamma(\frac{x}{beta}, gamm)} 
#' @usage precintcon.H(beta, gamm, pzero, x) 
#' @param beta is the \deqn{\beta} calculated by \deqn{\frac{\bar{x}}{\alpha}}.
#' @param gamm is the \deqn{\gamma} calculated by the \code{\link{precintcon.gamma}} function.
#' @param pzero is the probability of occur zeros in the serie.
#' @param x is the precipitation value.
#' @return The cumulative probability \code{H} defined by 
#' \deqn{pzero + (1 - pzero) \times pgamma(\frac{x}{beta}, gamm)} 
#' @seealso \code{\link{precintcon.ci.analysis}} 
#' @keywords precipitation concentration index 
precintcon.H <- function(beta, gamm, pzero, x) {
	if (x > 0)
		return(pzero + (1 - pzero) * pgamma(x/beta, gamm))
	else
		return(pzero)
}