#' guessdk_lik
#' @description Likelihood function for data with Don't Know. Used Internally.
#' @keywords internal
#' 
#' @param x    lgg, lgk, lgc, lkk, lcg, lck, and lcc
#' @param g1   guess
#' @param data transition matrix

guessdk_lik <- function(x, g1=x[8], data) 
{
	lgg <- x[1]
	lgk <- x[2] 
	lgc <- x[3]
	lkk <- x[4]
	lcg <- x[5]
	lck <- x[6]
	lcc <- x[7]
			
	vec <- NA
	vec[1] <- (1 - g1)*(1 - g1)*lgg
	vec[2] <- (1 - g1)*g1*lgg + (1 - g1)*lgk
	vec[3] <- (1 - g1)*lgc
	vec[4] <- (1 - g1)*g1*lgg
	vec[5] <- g1*g1*lgg+g1*lgk+lkk
	vec[6] <- g1*lgc
	vec[7] <- (1 - g1)*lcg
	vec[8] <- g1*lcg + lck
	vec[9] <- lcc
	
	-sum(data*log(vec))
}