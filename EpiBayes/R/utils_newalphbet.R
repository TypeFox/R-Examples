#' @title 
#' Computes Beta Parameters given a Mean and Variance
#' 
#' @description
#' This function is designed to approximate a density of a bounded (0, 1) function with a 
#'     beta distribution by equating the mean and variance to the first and second
#'     shape parameters. It does do by using the fact that:
#'     \deqn{\mu = \frac{\alpha}{\alpha + \beta}}{\mu = \alpha / (\alpha + \beta)}
#'     \deqn{\sigma^2 = \frac{\alpha \beta}{(\alpha + \beta)^2 + (\alpha + \beta + 1)}}{\sigma^2 = (\alpha \beta) / [(\alpha + \beta)^2 + (\alpha + \beta + 1)]}
#'
#' @param mu The mean of the function to be approximated.
#' @param sigma2 The variance of the function to be approximated.
#'
#' @return
#' The returned values are the first and second shape parameters of the fitted beta 
#'     distribution, stored in a 2 x 1 vector.
#'
#' @details
#' Some nonsensical answers could be returned if one simply uses the formulas for the 
#'     first and second shape parameters ourtight, especially concerning beta distributions
#'     at the extreme ends of the (0, 1) interval. To account for this, the function 
#'     handles two cases: (1) if the mean is greater than or equal to 1 or is less than or
#'     equal to 0, it is arbitrarily set to 0.99 or 0.01, respectively,  and (2)
#'     if either of the two shape parameters would be returned as a negative value, the 
#'     variance is increased for 1,000 iterations to try to force positive shape parameters.
#'     After those 1,000 iterations, an error is returned.
#'
#' @examples
#' ## Say we have that the tau posterior distribution from EpiBayes_ns() has mean and 
#' ## variance 0.01, and 0.015, respectively. The corresponding beta parameters will be:
#' utils_newalphbet(0.01, 0.015)
#'
#' ## If we provide a mean of 1, gives meaningful results
#' utils_newalphbet(1, 1)
#' 
#' ## If we provide a mean of 0, gives meaningful results
#' utils_newalphbet(0, 1)
#'
#' \dontrun{
#' ## Returns an error message if shape parameters cannot be coerced to be positive
#' utils_newalphbet(1, -1)
#' }
#' 
#' @export
#' @name utils_newalphbet
utils_newalphbet = function(mu, sigma2){
	
	## Set a mean of 1 or 0 to a mean of 0.99 or 0.01 to avoid nonsensical answers
	if (mu >= 1){
		mu = 0.99
	}
	
	if (mu <= 0){
		mu = 0.01
	}
		
	## Calculate the new shape parameters preliminarily	
	newalph.out = (mu^2 * (1 - mu) - sigma2 * mu) / sigma2
	newbet.out = (mu * (1 - mu)^2 - sigma2 * (1 - mu)) / sigma2
	
	## Set a counter for the while loop to avoid infinite loops
	i = 0
	
	## While loop cycles through and increases variance until both shape parameters are positive
	while (any(c(newalph.out, newbet.out) <= 0)){
		
		# Increase the variance by 10% each loop
		sigma2 = sigma2 * 0.1
		
		# Recalculate the shape parameters
		newalph.out = (mu^2 * (1 - mu) - sigma2 * mu) / sigma2
		newbet.out = (mu * (1 - mu)^2 - sigma2 * (1 - mu)) / sigma2
		
		# Error handling for negative shape parameters
		i = i + 1
		if (i > 1000) stop("TOO MANY ITERATIONS: the shape parameters remained negative with the given mean and variance after increasing the variance by 10% 1,000 times. Please provide a new mean and variance to yield positive shape parameters.")

	}
	
	## Calculate and return the new shape parameters
	newalph.out = (mu^2 * (1 - mu) - sigma2 * mu) / sigma2
	newbet.out = (mu * (1 - mu)^2 - sigma2 * (1 - mu)) / sigma2
	newalphbet.out = c(newalph.out, newbet.out)

	return(newalphbet.out)
}