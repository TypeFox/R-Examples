#' Hall-Buckley-Eagleson method
#'
#' Computes the cdf of a positively-weighted sum of chi-squared random variables with the Hall-Buckley-Eagleson (HBE) method.
#' @inheritParams sw
#' @keywords distribution
#' @references
#' \itemize{
#'   \item P. Hall. Chi squared approximations to the distribution of a sum of independent random variables. \emph{The Annals of Probability}, 11(4):1028-1036, 1983.
#'   \item M. J. Buckley and G. K. Eagleson. An approximation to the distribution of quadratic forms in normal random variables. \emph{Australian Journal of Statistics}, 30(1):150-159, 1988.
#' }
#' @export
#' @examples
#' #Examples taken from Table 18.6 in N. L. Johnson, S. Kotz, N. Balakrishnan. 
#' #Continuous Univariate Distributions, Volume 1, John Wiley & Sons, 1994.
#'
#' hbe(c(1.5, 1.5, 0.5, 0.5), 10.203)            # should give value close to 0.95
#' hbe(coeff=c(1.5, 1.5, 0.5, 0.5), x=10.203)    # specifying parameters
#' hbe(c(1.5, 1.5, 0.5, 0.5), c(0.627, 10.203))  # x is a vector, output approx. 0.05, 0.95

hbe <- function(coeff, x){
    if ( (missing(x)) || (missing(coeff)) )
        stop("missing an argument - need to specify \"coeff\" and \"x\"")

    if (checkCoeffsArePositiveError(coeff))
        stop(getCoeffError(coeff))

    if (checkXvaluesArePositiveError(x))
        stop(getXvaluesError(x))
    #end of error checking

	#compute cumulants and nu
	kappa <- c(sum(coeff), 2*sum(coeff^2), 8*sum(coeff^3) )
    K_1 <- sum(coeff)
    K_2 <- 2 * sum(coeff^2)
    K_3 <- 8 * sum(coeff^3)
	nu <- 8 * (K_2^3) / (K_3^2)
	
	#gamma parameters for chi-square
	gamma_k <- nu/2
	gamma_theta <- 2

	#need to transform the actual x value to x_chisqnu ~ chi^2(nu)
	#This transformation is used to match the first three moments
	#First x is normalised and then scaled to be x_chisqnu
	x_chisqnu_vec <- sqrt(2 * nu / K_2) * (x - K_1) + nu
	
	#now this is a chi_sq(nu) variable
	p_chisqnu_vec <- pgamma(x_chisqnu_vec, shape=gamma_k, scale=gamma_theta)
	return(p_chisqnu_vec)
}

