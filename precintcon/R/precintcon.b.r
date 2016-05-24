#' @noRd
#' @name precintcon.b 
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.b 
#' @title b 
#' @description Calculates b via the least-squares method. 
#' @usage precintcon.b(X, Y) 
#' @param X is the cumulative percentage of rainy days.
#' @param Y is the cumulative percentage of rainfall amounts.
#' @return \code{b} 
#' @seealso \code{\link{ci}} 
#' @keywords precipitation concentration index 
precintcon.b <- function(X, Y) {
	
    N <- length(X)
	
	if (length(X) == length(Y))
		return(
			(
			   (N * sum(X * log(Y, base=exp(1)))) + 
			   (sum(X) * sum(log(X,base=exp(1)))) - 
			   (N * sum(X * log(X,base=exp(1)))) - 
			   (sum(X) * sum(log(Y,base=exp(1))))
			) / 
			(N * sum(X^2) - (sum(X)^2))/100)
	
	else
		stop("X and Y have differents lengths!")
}
