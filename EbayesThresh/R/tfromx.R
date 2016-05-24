"tfromx" <-
function(x, prior = "laplace", bayesfac = FALSE, a = 0.5)
{
#  given the data x, the prior, and any other parameters, 
#   find the threshold
#   corresponding to the marginal maximum likelihood estimator
#   of the mixing weight.
#
	if ( prior=="laplace" && is.na (a) ) 
	{ wa  <-  wandafromx( x)
	w  <-  wa$w
	a  <-  wa$a 
	}  	else	
	{w <- wfromx(x, prior, a = a)}
	t <- tfromw(w, prior, a = a, bayesfac = bayesfac)
	return(t)
}
