lognormalEB<-function(Observed, Expected, maxiter = 20, tol = 1e-05)
{
	n<-length(Observed)

	#Initial values
	b<-log((Observed+.5)/Expected)
	m0<-mean(b)
	v0<-var(b)

	m1<-m0
	v1<-(v0*sum( 1/(1+v0*(Observed+.5)) ) +sum( (b-m1)*(b-m1) ) )/n

	b<-m1+(Observed+.5)*v1*log((Observed+.5)/Expected)-v1/2
	b<-b/(1+(Observed+.5)*v1)

	iter<-1
	while( ((abs((m0-m1)/(m0+m1))>tol) || (abs((v0-v1)/(v0+v1))>tol)) &&(iter<=maxiter))
	{
		m0<-m1
		v0<-v1

		m1<-mean(b)
		v1<-(v0*sum(1/(1+v0*(Observed+.5))) +sum((b-m1)*(b-m1)))/n

		b<-m1+(Observed+.5)*v1*log((Observed+.5)/Expected)-v1/2
		b<-b/(1+(Observed+.5)*v1)

		iter<-iter+1
	}

#	print(iter)
	b

	return(list(n = length(Observed), phi = m1, sigma2=v1,
        smthrr = b))
}
