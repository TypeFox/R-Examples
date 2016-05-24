#' @noRd 
#' @name precintcon.ln.a
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.ln.a 
#' @title ln(a)
#' @description Performes the calculation of ln(a) via the least-squares method. 
#' @usage precintcon.ln.a(X, Y) 
#' @param X is the cumulative percentage of rainy days.
#' @param Y is the cumulative percentage of rainfall amounts.
#' @return ln(a) 
#' @seealso \code{\link{precintcon.ci.analysis}} 
#' @keywords precipitation concentration index 
precintcon.ln.a <- function(X, Y) {
	
	N <- length(X)
	
	if (length(X) == length(Y))
       return(
	      (
		     (sum(X^2) * sum(log(Y,base=exp(1)))) + 
		     (sum(X) * sum(X*log(X, base=exp(1)))) - 
		     (sum(X^2) * sum(log(X,base=exp(1)))) - 
		     (sum(X) * sum(X * log(Y,base=exp(1))))
		  ) / 
	      (N * sum(X^2) - (sum(X)^2))
       )
   
   else
	   writeLines("X and Y with differents lengths!")
}

