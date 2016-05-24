"beta.cauchy" <-
function(x)
{
#
#   Find the function beta
#    for the mixed normal prior with Cauchy tails.  It is assumed that the 
#    noise variance is equal to one.  
#
	phix <- dnorm(x)
	j <- (x != 0)
	beta <- x
	beta[!j] <- -1/2
	beta[j] <- (dnorm(0)/phix[j] - 1)/x[j]^2 - 1
	return(beta)
}
