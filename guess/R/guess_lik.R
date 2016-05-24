#' guess_lik
#' @description Likelihood function for data without Don't Know. Used Internally.
#' @keywords internal
#' 
#' @param x    lgg, lgk, lkk
#' @param g1   guess
#' @param data transition matrix

guess_lik <- function(x, g1=x[4], data) 
{
	lgg <- x[1]
 	lgk <- x[2]
 	lkk <- x[3]
			
	vec <- NA
 	vec[1] <- (1-g1)*(1-g1)*lgg
 	vec[2] <- (1-g1)*g1*lgg + (1-g1)*lgk
 	vec[3] <- (1-g1)*g1*lgg
 	vec[4] <- g1*g1*lgg+g1*lgk+lkk
 	
 	-sum(data*log(vec))	
}