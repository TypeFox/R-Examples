kn.iscluster<-function(data, idx, idxorder, alpha, fractpop, use.poisson=TRUE, 
	model="poisson", R, mle, ...)
{
        #Indexes of regions in the balls ordered by distance to the centre
        localidx<-idxorder[idx[idxorder]]


	knboot<-switch(model,
	permutation = boot(data[idxorder,], statistic=kullnagar.boot, R=R, fractpop=fractpop, use.poisson=use.poisson, ...),
	multinomial=boot(data[localidx,], statistic=kullnagar.pboot, sim="parametric", ran.gen=multinom.sim,  R=R, fractpop=fractpop, use.poisson=use.poisson, mle=list(n=mle$n, p=mle$p[localidx]), ... ),
	poisson = boot(data[localidx,], statistic=kullnagar.pboot, sim="parametric", ran.gen=poisson.sim,  R=R, fractpop=fractpop, use.poisson=use.poisson, mle=list(n=mle$n, lambda=mle$lambda[localidx]), ... ),
	negbin = boot(data[localidx,], statistic=kullnagar.pboot, sim="parametric", ran.gen=negbin.sim,  R=R, fractpop=fractpop, use.poisson=use.poisson, mle=list(n=mle$n, size=mle$size,prob=mle$prob[localidx]), ... )
	)

	if(is.null(knboot$t0))
		return(c(NA, NA, NA, NA))

	pvalue<-(sum(knboot$t[,1]>knboot$t0[1])+1)/(R+1)
	return(c( knboot$t0[1], alpha>pvalue, pvalue, knboot$t0[2]) )
}


#Compute the statistic around the first location(row) in the dataframe
#Data must be order according to distance to the center of the circle
kullnagar.stat<-function(data, fractpop, use.poisson=TRUE, log.v=FALSE)
{

	n<-length(data[[1]])

	if(use.poisson)
	{
	#	r<-kullnagar.stat.poisson(data, fractpop, log.v=FALSE)
		r<-.Call("Rkn_poisson", as.numeric(data$Observed),
		data$Expected, fractpop, PACKAGE="DCluster")

		if(!log.v)
			r[1]<-exp(r[1])

		r<-c(value=r[1], size=r[2])
	}
	else
		r<-kullnagar.stat.bern(data, fractpop, log.v)

	return(r)
}


#Bernouilli version
kullnagar.stat.bern<-function(data, fractpop, log.v=FALSE)
{
	n<-length(data[[1]])


	csumpop<-cumsum(data$Population)
	P<-csumpop[n]

	p<-sum(csumpop<(fractpop*P))+1
	if(p>=n) p<-n-1 #The ratio is 1 when all regions are considered

	L<-0
	O<-sum(data$Observed)
	csumobs<-cumsum(data$Observed[1:p])

	#log(L_0)
	l0<-O*log(O)+(P-O)*log(P-O)-P*log(P)

	#Size of the cluster, number of regions from the centre
	size<-1
	for(i in 1:p)
	{
		if( (csumobs[i]*(P-csumpop[i])) > (csumpop[i]*(O-csumobs[i])) )
		{
			#Likelihood inside ball Z
			l<-csumobs[i]*log(csumobs[i])
			
			aux<-csumpop[i]-csumobs[i]
			l<-l+ aux*log(aux) 
			
			l<-l- csumpop[i]*log(csumpop[i])
				
			
			#Likelihood outside ball Z
			aux<-O-csumobs[i]
			l<-l+aux*log(aux)
			
			aux<-P-csumpop[i]-O+csumobs[i]
			l<-l+aux*log(aux)
			
			aux<-P-csumpop[i]
			l<-l-aux*log(aux)

			#L_0
			l<-l-l0

		}
		else
		{
			l<-0
		}

		if(l>L)
		{
			L<-l
			size<-i
		}
	}

	if(!log.v) L<-exp(L)

	return(c(value=L, size=size))
}

#Poisson version
kullnagar.stat.poisson<-function(data, fractpop, log.v=FALSE)
{
	n<-length(data[[1]])

	csumexp<-cumsum(data$Expected)
	E<-csumexp[n]

	p<-sum(csumexp<(fractpop*E))+1
	if(p>=n) p<-n-1 #The ratio is 1 when all regions are considered

	L<-0
	O<-sum(data$Observed)
	csumobs<-cumsum(data$Observed[1:p])


	#log(L_0)
	l0<-O*log(O/E)

	#Size of the cluster, number of regions from the centre
	size<-1
	for(i in 1:p)
	{
		difobs<-O-csumobs[i]
		difexp<-E-csumexp[i]

		if( (csumobs[i]*difexp) > (csumexp[i]*difobs) )
		{
			#Likelihood inside ball Z
			l<-csumobs[i]*log(csumobs[i]/csumexp[i])
			
			#Likelihood outside ball Z
			l<-l+difobs*log(difobs/difexp)

			#Divide by L_0
			l<-l-l0
		}
		else
		{
			l<-0
		}

		if(l>L)
		{
			L<-l
			size<-i
		}
	}

	if(!log.v) L<-exp(L)

	return(c(value=L, size=size))
}
