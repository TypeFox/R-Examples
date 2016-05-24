pmle.exp <-
function(x,m0=1,lambda=0,inival=NULL,len=10,niter=50,tol=1e-6,rformat=FALSE)
#This function compute the PMLE of parameters under a exponential mixture
#x: 		data, can be either a vector or a matrix with the 1st column being the observed values 
#        		and the 2nd column being the corresponding frequencies. 
#m0: 		order of finite mixture model.
#inival:	initial values for the EM-algorithm
#lambda:	size of penalty function of mixing proportions.
#len: 	number of initial values chosen for the EM-algorithm.
#niter:     least number of iterations for all initial values in the EM-algorithm.
#tol: 	tolerance value for the convergence of the EM-algorithm.
#rformat	format for output, rformat=T means the format of output is determined by R software.
#		rformat=F means the format of output is determined by our default setting. When the output is
#		larger than 0.001, it is determined by round(output,3); When the output is less than 0.001,
#		it is determined by signif(output,3).
{
	if (is.data.frame(x))
	{	
		if (ncol(x)==2)
			x=as.matrix(x)
		if (ncol(x)==1 | ncol(x)>2)
			x=x[,1]
	}
	if (is.matrix(x))
	{
		xx=c()
		for (i in 1:nrow(x))
			xx=c(xx,rep(x[i,1],x[i,2]))
		x=xx
	}

	out=phi0.exp(x,m0,lambda,inival,len,niter,tol)

	alpha=out$alpha
	theta=out$theta
	loglik=out$loglik
	ploglik=out$ploglik

	if (rformat==F)
	{
		alpha=rousignif(alpha)
		theta=rousignif(theta)
		loglik=rousignif(loglik)
		ploglik=rousignif(ploglik)
	}

	if (lambda==0 | m0==1)
		list('MLE of mixing proportions:'=alpha,
		'MLE of component parameters:'=theta,
		'log-likelihood:'=loglik)
	else
		list('PMLE of mixing proportions:'=alpha,
		'PMLE of component parameters:'=theta,
		'log-likelihood:'=loglik,
		'Penalized log-likelihood:'=ploglik)
}
