phi.norm <-
function(x,m0,lambda,inival,len,niter,tol)
{
#x: 		data, can be either a vector or a matrix with the 1st column being the observed values
#       		and the 2nd column being the corresponding frequencies. 
#m0:	 	order of finite normal mixture model.
#lambda:	size of penalized function of mixing distribution
#inival:	initial values chosen for the EM-algorithm
#len:		number of initial values chosen for the EM-algorithm.	
#niter:     least number of iterations for all initial values in the EM-algorithm.
#tol:		tolerance value for the convergence of the EM-algorithm.
	if (m0>1)
	{
		sn=var(x)
		dn=sd(x)
		maxx=max(x)
		minx=min(x)
		n=length(x)
		
		if (is.vector(inival))
			inival=t(inival)
		if (is.null(inival)==F)
			len=nrow(inival)
		output=c()
		for (i in 1:len)
		{
			###initial values for EM-algorithm
			if (is.null(inival))
			{
				###random initial values
				alpha=runif(m0,0,1)
				alpha=alpha/sum(alpha)
				mu=runif(m0,minx,maxx)	
				sigma=runif(m0,0.25,2)*sn
			}
			else
			{	
				###given initial values in input
				alpha=inival[1:m0]
				mu=inival[(m0+1):(2*m0)]
				sigma=sqrt(inival[(2*m0+1):(3*m0)])
			}
			para0=c(alpha,mu,sigma)

			for (j in 1:niter)###run niter EM-iterations first
			{
				outpara=iter1.norm(x,para0,lambda)
				para0=outpara[1:(3*m0)]
			}
			output=rbind(output,outpara)	
		}
		index=which.max(output[,(3*m0+2)])
		para0=output[index,1:(3*m0)]
		pln0=output[index,(3*m0+2)]
		err=1
		t=0
		while (err>tol & t<2000)###EM-iteration with the initial value with the largest penalized log-likelihood
		{
			outpara=iter1.norm(x,para0,lambda)
			para0=outpara[1:(3*m0)]
			pln1=outpara[3*m0+2]
			err=pln1-pln0
			pln0=pln1
			t=t+1
		}
		list('alpha'=outpara[1:m0],
		'means'=outpara[(m0+1):(2*m0)],
		'variances'=outpara[(2*m0+1):(3*m0)],
		'loglik'=outpara[3*m0+1],
		'ploglik'=outpara[3*m0+2])
	}
	else 
	{
		mu0=mean(x)
		sig0=var(x)
		ln=sum(log(dnorm(x,mu0,sqrt(sig0))))
		list('alpha'=1,'means'=mu0,'variances'=sig0,'loglik:'=ln,'ploglik:'=ln)
	}
}
