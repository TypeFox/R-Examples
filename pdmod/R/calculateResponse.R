##' Calculate response from the estimate
##' 
##' Given an estiamtes probability of reward beween 0 and 1, calculates a response rate 
##' (i.e. the measured response of the animal such as visits to the food delivery system)
##' 
##' @param k		Response rate parameter
##' @param rmax		Maximum response
##' @param est		Vector of estimates 
##' @return 		Vector of responses 
##' @export
##' @seealso \code{\link{Constants}}, \code{\link{isTimedVector}}, \code{\link{verifyTimedVector}}
##' @author Chloe Bracis
##' @examples
##' calculateResponse(0.8, 10, runif(20))

calculateResponse = function( k, rmax, est ) rmax * est / ( est + k * ( 1 - est ) )
