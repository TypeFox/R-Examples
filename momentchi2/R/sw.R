#' Satterthwaite-Welch method
#'
#' Computes the cdf of a positively-weighted sum of chi-squared random variables with the Satterthwaite-Welch (SW) method.
#' @param coeff The coefficient vector. All values must be greater than 0.
#' @param x The vector of quantile values. All values must be greater than 0.
#' @keywords distribution
#' @references 
#' \itemize{ 
#'  \item B. L.Welch. The significance of the difference between two means when the population variances are unequal. \emph{Biometrika}, 29(3/4):350-362, 1938. 
#'  \item F. E. Satterthwaite. An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2(6):110-114, 1946.
#'  \item G. E. P. Box Some theorems on quadratic forms applied in the study of analysis of variance problems, I. Effects of inequality of variance in the one-way classification. \emph{The Annals of Mathematical Statistics}, 25(2):290-302, 1954.
#' }    
#' @export
#' @examples
#' #Examples taken from Table 18.6 in N. L. Johnson, S. Kotz, N. Balakrishnan. 
#' #Continuous Univariate Distributions, Volume 1, John Wiley & Sons, 1994.
#'
#' sw(c(1.5, 1.5, 0.5, 0.5), 10.203)            # should give value close to 0.95
#' sw(coeff=c(1.5, 1.5, 0.5, 0.5), x=10.203)    # specifying parameters
#' sw(c(1.5, 1.5, 0.5, 0.5), c(0.627, 10.203))  # x is a vector, output close to 0.05, 0.95
sw <- function(coeff, x){

    if ( (missing(x)) || (missing(coeff)) )
        stop("missing an argument - need to specify \"coeff\" and \"x\"")

    if (checkCoeffsArePositiveError(coeff))
        stop(getCoeffError(coeff))

    if (checkXvaluesArePositiveError(x))
        stop(getXvaluesError(x))
    #end of error checking

#checkCoeffsArePositive(coeff)
#   checkXvaluesArePositive(x)
# try(checkCoeffsArePositive(coeff), silent=FALSE)
#try(checkXvaluesArePositive(x), silent=FALSE)
	
	#compute cumulant and ratio of cumulants
	w_val <- sum(coeff)
	u_val <- sum(coeff^2) / (w_val^2)
	
	#now the G k and theta:
	gamma_k <- 0.5 / u_val
	gamma_theta <- 2 * u_val*w_val
	
	#the actual x value
	p_sw <- pgamma(x, shape=gamma_k, scale=gamma_theta)	
	return(p_sw)
}

