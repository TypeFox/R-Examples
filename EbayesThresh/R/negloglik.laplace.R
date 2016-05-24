"negloglik.laplace" <-
function(xpar, xx)
{
#
#  marginal negative log likelihood function for laplace prior
#  
#  xx data
#  xpar vector of two parameters:
#      xpar[1] :  threshold
#      xpar[2] :  scale factor a
#
	a <- xpar[2]
	w <- wfromt(xpar[1], a = a)
	loglik <- sum(log(1 + w * beta.laplace(xx, a)))
	return( - loglik)
}
