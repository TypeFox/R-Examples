#' @title \code{dkay} The density function of the K distribution
#'
#' @description The K density function on \code{df} degrees of freedom and non-centrality parameter \code{ncp}.  
#' 
#'    A K distribution is the square root of a chi-square divided by its degrees of freedom.  That is, if x is chi-squared on m degrees of freedom, then y = sqrt(x/m) is K on m degrees of freedom.
#'    Under standard normal theory, K is the distribution of the pivotal quantity  s/sigma where s is the sample standard deviation and sigma is the standard deviation parameter of the normal density.  K is the natural distribution for tests and confidence intervals about sigma.
#'    K densities are more nearly symmetric than are chi-squared and concentrate near 1.  As the degrees of freedom increase, they become more symmetric, more concentrated, and more nearly normally distributed. 
#'     
#'     
#'
#' @export dkay
#' 
#' @param x A vector of values at which to calculate the density.
#' @param df Degrees of freedom (non-negative, but can be non-integer).
#' @param ncp Non-centrality parameter (non-negative).
#' @param log.p logical; if \code{TRUE}, probabilities are given as log(p).
#'
#' @return \code{dkay} gives the density evaluated at the values of \code{x}.
#' 
#' Invalid arguments will result in return value NaN, with a warning. 
#' 
#' The length of the result is the maximum of the lengths of the numerical arguments for the other functions.
#' 
#' The numerical arguments are recycled to the length of the result. Only the first elements of the logical arguments are used.
#' 
#'
#' @note All calls depend on analogous calls to chi-squared functions.  See \code{dchisq} for details on non-centrality parameter calculations.
#'
#'
#' @examples
#' 
#' dkay(1, 20)
#' #
#' # compare K density to that of chi as degrees of freedom increase
#' op <-par(mfrow=c(1,2))
#' p <- seq(0.001, .999, 0.001)
#' #
#' # First get all the chi-square densities and plot them
#' xchi5 <- qchisq(p,5)
#' dchi5 <- dchisq(xchi5,5)
#' xchi10 <- qchisq(p,10)
#' dchi10 <- dchisq(xchi10,10)
#' xchi20 <- qchisq(p,20)
#' dchi20 <- dchisq(xchi20,20)
#' xchi30 <- qchisq(p,30)
#' dchi30 <- dchisq(xchi20,30)
#' xlim <- range(xchi5, xchi10, xchi20, xchi30)
#' ylim <- range(dchi5, dchi10, dchi20, dchi30)
#' plot(xchi5, dchi5, type="l", xlab="x", ylab="density", 
#'      xlim=xlim, ylim=ylim,  
#'      main="chi-squared densities")
#' lines(xchi10, dchi10, lty=2)
#' lines(xchi20, dchi20, lty=3)
#' lines(xchi20, dchi30, lty=4)
#' legend("topright",  
#'        legend=c("df = 5", "df = 10", "df = 20", "df = 30"),  
#'        lty=c(1,2,3,4),  
#'        title="degrees of freedom",  
#'        cex=0.75, bty="n")
#' #
#' # Now get all the K densities and plot them
#' xkay5 <- qkay(p,5)
#' dkay5 <- dkay(xkay5,5)
#' xkay10 <- qkay(p,10)
#' dkay10 <- dkay(xkay10,10)
#' xkay20 <- qkay(p,20)
#' dkay20 <- dkay(xkay20,20)
#' xkay30 <- qkay(p,30)
#' dkay30 <- dkay(xkay20,30)
#' xlim <- range(xkay5, xkay10, xkay20, xkay30)
#' ylim <- range(dkay5, dkay10, dkay20, dkay30)
#' plot(xkay5, dkay5, type="l",  
#'      xlab="x", ylab="density", 
#'      xlim=xlim, ylim=ylim,  
#'      main="K densities")
#' lines(xkay10, dkay10, lty=2)
#' lines(xkay20, dkay20, lty=3)
#' lines(xkay20, dkay30, lty=4)
#' legend("topright",  
#'        legend=c("df = 5", "df = 10", "df = 20", "df = 30"),  
#'        lty=c(1,2,3,4),  
#'        title="degrees of freedom",  
#'        cex=0.75, bty="n")
#' par(op)
#' #
#' # Note how K densities are more nearly symmetric and concentrate near 1.  
#' # As the degrees of freedom increase, they become more symmetric,  
#' # more concentrated, and more nearly normally distributed.
dkay <- function(x, df, ncp=0, log.p = FALSE) {
				chi <- df * x^2
				chincp <- df * ncp^2
				fchi <- dchisq(chi, df, chincp, log.p)
				2 * df * x * fchi
				} 