emtest.norm <-
function(x,m0=1,pens=NULL,inival=NULL,len=10,niter=50,tol=1e-6,k=3,rformat=FALSE)
{
#x: 		data, can be either a vector or a matrix with the 1st column being the observed values
#       		and the 2nd column being the corresponding frequencies. 
#m0:	 	order of finite normal mixture model under the null hypothesis.
#pens:	a 2-dimensions vector being the size of penalized functions for mixing propotion and variance.
#inival:	initial values chosen for the EM-algorithm to compute the PMLE under the null model
#len:		number of initial values chosen for the EM-algorithm.
#niter:     least number of iterations for all initial values in the EM-algorithm.
#tol:		tolerance value for the convergence of the EM-algorithm.
#k:		number of EM iterations to obtain EM-test statistic.	
#rformat	format for output, rformat=T means the format of output is determined by R software.
#		rformat=F means the format of output is determined by our default setting. When the output is
#		larger than 0.001, it is determined by round(output,3); When the output is less than 0.001,
#		it is determined by signif(output,3).
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

	###PMLE and log-likelihood value under the null model	
	para0=phi.norm(x,m0,0,inival,len,niter,tol)
	out.null=para0[[4]]

	###Size of penalized function
	if (is.null(pens))
	{
		pens[1]=1
		if (m0==1)
			pens[2]=0.25
		if (m0==2)
			pens[2]=formula2(para0,length(x))
		if (m0==3)
			pens[2]=formula3(para0,length(x))
		if (m0>3)
			pens[2]=0.2
	}

	###create beta vectors for EM-test
	bbeta=c()
	for(h in 1:m0) 
		bbeta=rbind(cbind(bbeta,rep(0.1,3^{h-1})),
            	cbind(bbeta,rep(0.3,3^{h-1})),
            	cbind(bbeta,rep(0.5,3^{h-1})))
	output=c()
	for (i in 1:(3^m0))
	{
		beta=bbeta[i,]
		###PMLE given beta
		para1=maxmm.norm(x,beta,pens[2],para0,len,niter,tol)
		###k iterations of EM-algorithm
		out=emiter.norm(x,para1,beta,pens,para0,k)
		output=rbind(output,out)
	}
	index=which.max(output[,(6*m0+1)])
	outpara=output[index,1:(6*m0)]
	out.alt=output[index,(6*m0+1)]

	###EM-test statistic
	emstat=2*(out.alt-out.null)
	emstat=as.matrix(emstat)
	rownames(emstat)=NULL
	emstat=as.vector(emstat)
	
	###output
	alpha=para0[[1]]
	mean=para0[[2]]
	variance=para0[[3]]
	t0=rbind(alpha,mean,variance)
	alpha=outpara[1:(2*m0)]
	mean=outpara[(2*m0+1):(4*m0)]
	variance=outpara[(4*m0+1):(6*m0)]
	t1=rbind(alpha,mean,variance)
	p=pchisq(emstat,(2*m0),lower.tail=F)	

	if (rformat==F)
	{
		t0=rousignif(t0)
		t1=rousignif(t1)
		emstat=rousignif(emstat)
		p=rousignif(p)		
		pens=rousignif(pens)	
	}
	
	list(
	'MLE of Parameters under null hypothesis (order = m0)'=t0,
	'Parameter estimates under the order = 2m0'=t1,
	'EM-test Statistics'=emstat,
	'P-values'=p,
	'Level of penalty'=pens)
}
