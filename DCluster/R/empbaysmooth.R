empbaysmooth<-function(Observed, Expected, maxiter=20, tol=1e-5)
{
	if(length(Observed)!=length(Expected))
	{
		print("Lengths of the two vectors differs")
		return(NULL)
	}
	
	n<-length(Observed)
	idx<-Expected>0
	n1<-length(Expected[idx])#Number of non-zero values
	
	#Starting point: the smoothed R.R. are Observed/Expected
	smthrr<-rep(0,n)
	smthrr[idx]<-Observed[idx]/Expected[idx]
	
	m0<-mean(smthrr[idx])
	v0<-var(smthrr[idx])

	if(v0==0)#No variability, i.e., Observed = K
	{
		print("Observed cases are equal to a constant.")
		return( list(nu=NA, alpha=NA, smthrr=rep(NA,n), niter=0) )
	}

	#Initial values for the gamma parameters
	nu<-m0*m0/v0
	alpha<-m0/v0

	smthrr[idx]<-(nu+Observed[idx])/(alpha+Expected[idx])

	m<-mean(smthrr[idx])
	v<- sum( (1+alpha/Expected[idx]) * ((smthrr[idx]-m)^2) )/(n1-1)

	iter<-1
	while(  (  ( abs(m-m0) >tol*(m+m0) ) || ( abs(v-v0) >tol*(v+v0) ) )   && ( iter<=maxiter) )
	{
		#Updated values for the gamma parameters
		nu<-m*m/v
		alpha<-m/v

		smthrr[idx]<-(Observed[idx]+nu)/(Expected[idx]+alpha)
		
		#Previous mean and variance
		m0<-m
		v0<-v

		m<-mean(smthrr[idx])
		v<- sum( (1+alpha/Expected[idx]) * ((smthrr[idx]-m)^2) )/(n1-1)

		iter<-iter+1
	}


	return( list(n=length(Observed), nu=nu, alpha=alpha, smthrr=smthrr, niter=iter) )

}
