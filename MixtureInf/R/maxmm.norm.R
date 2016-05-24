maxmm.norm <-
function(x,beta,an,para0,len,niter,tol)
{
	dn=sd(x)
	m=length(beta)
	alpha0=para0[[1]]
	mu0=para0[[2]]
	sigma0=para0[[3]]
	output=c()

	###Calculate the cut points of parameter space of mu
	eta=rep(0,(m+1))
	eta[1]=min(x)
	eta[m+1]=max(x)
	if (m>1)
		eta[2:m]=(mu0[1:(m-1)]+mu0[2:m])/2
	
	for (i in 1:len)
	{
		###initial values for EM-algorithm
		alpha=c()
		mu=c()
		sigma=c()
		for (j in 1:m)
		{	
			alpha=c(alpha,alpha0[j]*beta[j],alpha0[j]*(1-beta[j]))
			mu=c(mu,runif(2,eta[j],eta[j+1]))
			sigma=c(sigma,runif(2,0.25*sigma0[j],2*sigma0[j]))
		}
		para=c(alpha,mu,sigma)

		for (j in 1:niter)###run niter EM-iterations first
		{
			outpara=iter2.norm(x,para,beta,an,para0)
			para=outpara[1:(6*m)]
		}
		output=rbind(output,outpara)	
	}
	index=which.max(output[,(6*m+1)])
	para=output[index,1:(6*m)]
	pln0=output[index,(6*m+1)]
	err=1
	t=0
	while (err>tol & t<2000)###EM-iteration with the initial value with the largest penalized log-likelihood
	{
		outpara=iter2.norm(x,para,beta,an,para0)
		para=outpara[1:(6*m)]
		pln1=outpara[6*m+1]
		err=pln1-pln0
		pln0=pln1
		t=t+1
	}
	para
}
