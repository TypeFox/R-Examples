#' Weighted mean, variance and standard deviation calculations
#' 
#' \code{wt.mean} calculates the mean given a weighting of the values. \cr \cr
#' \code{wt.var} is the unbiased variance of the weighted mean calculation
#' using equations of GNU Scentific Library
#' (\url{http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.htmland}.\cr\cr
#' \code{wt.sd} is the standard deviation of the weighted mean calculated as
#' the sqrt of \code{wt.var}. \cr \cr \bold{Note:} NA data is automatically
#' ommitted from analysis.
#' 
#' 
#' @param x is a vector of numerical data.
#' @param wt is a vector of equal length to \code{x} representing the weights.)
#' @return returns a single value from analysis requested.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' #define simple data
#' x = 1:25 # set of numbers
#' wt = runif(25) #some arbitrary weights
#' 
#' #display means & variances (unweighted and then weighted)
#' mean(x); wt.mean(x,wt)
#' var(x); wt.var(x,wt)
#' sd(x); wt.sd(x,wt)
#' 
#' 
#' @export 
wt.mean <- function(x,wt) {
	s = which(is.finite(x*wt)); wt = wt[s]; x = x[s] #remove NA info
	return( sum(wt * x)/sum(wt) ) #return the mean
}

#' @rdname wt.mean
#' @export
wt.var <- function(x,wt) {
	s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
	xbar = wt.mean(x,wt) #get the weighted mean
	return( sum(wt *(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2))) ) #return the variance
} 

#' @rdname wt.mean
#' @export
wt.sd <- function(x,wt) { 
	return( sqrt(wt.var(x,wt)) ) #return the standard deviation
} 
