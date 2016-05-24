"postmean" <-
function(x, w, prior = "laplace", a = 0.5)
{
#
# find the posterior mean for the appropriate prior for 
#   given x and w and a, assuming the error variance
#   is 1.  
#
	pr <- substring(prior, 1, 1)
	if(pr == "l")
		mutilde <- postmean.laplace(x, w, a = a)
	if(pr == "c")
		mutilde <- postmean.cauchy(x, w)
	return(mutilde)
}
