"tfromw" <-
function(w, prior = "laplace", bayesfac = FALSE, a = 0.5)
{
#  given the vector of weights w, find the threshold or vector of
#   thresholds corresponding to these weights, under the specified prior.
#
#  if bayesfac=TRUE the Bayes factor thresholds are found, otherwise the posterior median
#   thresholds are found.
#
#  if the Laplace prior is used, a gives the value of the scale factor
#
	pr <- substring(prior, 1, 1)
	if(bayesfac) {
		z <- 1/w - 2
		if(pr == "l")
			tt <- vecbinsolv(z, beta.laplace, 0, 10, a = a)
		if(pr == "c")
			tt <- vecbinsolv(z, beta.cauchy, 0, 10)
	}
	else {
		zz <- rep(0, length(w))
		if(pr == "l")
			tt <- vecbinsolv(zz, laplace.threshzero, 0, 10, w = w, 
				a = a)
		if(pr == "c")
			tt <- vecbinsolv(zz, cauchy.threshzero, 0, 10, w = w)
	}
	return(tt)
}
