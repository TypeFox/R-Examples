#' Wood's F method
#'
#' Computes the cdf of a positively-weighted sum of chi-squared random variables with the Wood F (WF) method.
#' @inheritParams sw
#' @details Note that there are pathological cases where, for certain coefficient vectors (which result in certain cumulant values), the Wood F method will be unable to match moments (cumulants) with the three-parameter \emph{F} distribution. In this situation, the HBE method is used, and a warning is displayed. A simple example of such a pathological case is when the coefficient vector is of length 1. Note that these pathological cases are rare; see (Wood, 1989) in the references.
#' @keywords distribution
#' @references
#' \itemize{
#'  \item A. T. A. Wood. An F approximation to the distribution of a linear combination of chi-squared variables. \emph{Communications in Statistics-Simulation and Computation}, 18(4):1439-1456, 1989.
#' }
#' @export
#' @examples
#' #Examples taken from Table 18.6 in N. L. Johnson, S. Kotz, N. Balakrishnan. 
#' #Continuous Univariate Distributions, Volume 1, John Wiley & Sons, 1994.
#'
#' wf(c(1.5, 1.5, 0.5, 0.5), 10.203)            # should give value close to 0.95
#' wf(coeff=c(1.5, 1.5, 0.5, 0.5), x=10.203)    # specifying parameters
#' wf(c(1.5, 1.5, 0.5, 0.5), c(0.627, 10.203))  # x is a vector, output approx. 0.05, 0.95
#' wf(c(0.9), 1)                                # pathological case, warning, uses hbe()

wf <- function(coeff, x){
## warning: pathological case, uses hbe()
    if ( (missing(x)) || (missing(coeff)) )
        stop("missing an argument - need to specify \"coeff\" and \"x\"")

    if (checkCoeffsArePositiveError(coeff))
        stop(getCoeffError(coeff))

    if (checkXvaluesArePositiveError(x))
        stop(getXvaluesError(x))
    #end of error checking


	#r=1: 2^0*0! = 1
	K_1 <- sum(coeff)
	#r=2: 2^1*1! = 2
	K_2 <- 2*sum(coeff^2)
	#r=3: 2^2 *2! = 4*2=8
	K_3 <- 8*sum(coeff^3)
	
	#now need to calculate r's, in order to get beta
	r_1 <- 4 * K_2^2 * K_1  +  K_3 * (K_2 - K_1^2)
	r_2 <- K_3 * K_1 - 2* K_2^2
	
	if ((r_1==0) || (r_2==0)){
		warning("Either r1 or r2 equals 0: running hbe instead.")
		p_vec <- hbe(coeff, x)
	} else{
		#do the rest of wf
		beta <- r_1/r_2
		alpha_1 <- 2*K_1 * (K_3*K_1  +  K_1^2 * K_2 - K_2^2)/r_1
		alpha_2 <- 3 + 2*K_2*(K_2 + K_1^2)/r_2
		xvec_rescaled <- x * alpha_2/(alpha_1*beta)
		p_vec <- pf(xvec_rescaled, df1=2*alpha_1, df2=2*alpha_2)
	}
	return(p_vec)
}
