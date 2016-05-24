bn.iscluster<-function(data, idx, idxorder, alpha, k, model="poisson", R=999, mle)
{
        #Indexes of regions in the balls ordered by distance to the centre
	localidx<-idxorder[idx[idxorder]]

		
	bnboot<-switch(model,
	permutation = boot(data[localidx,], statistic=besagnewell.boot, R=R, k=k),
	multinomial = besagnewell.stat(data[localidx,], k ),
	poisson = besagnewell.stat(data[localidx,], k ),
	negbin = boot(data[localidx, ], statistic=besagnewell.pboot, R=R, sim="parametric",ran.gen=negbin.sim,  mle=list(n=sum(idx),size=mle$size,prob=mle$prob[localidx]), k=k)
	)

	stat<-ifelse(model=="permutation" | model=="negbin", bnboot$t0[1], bnboot[1])
	
	if(is.null(stat))
		return(c(NA,NA,NA,NA))

	pvalue<-switch(model,
	permutation=(1+.5*sum(bnboot$t[,1]==stat)+.5*sum(bnboot$t[,1]<stat))/(R+1),
	multinomial=1-pbinom(k-1, size=mle$n, prob=sum(mle$p[localidx[1:stat]]) ),
	poisson=1 - ppois(k-1, sum(mle$lambda[localidx[1:stat]])),
	negbin=(1+.5*sum(bnboot$t[,1]==stat)+.5*sum(bnboot$t[,1]<stat))/(R+1)
	)
		
	return(c(stat, (alpha>pvalue), pvalue, stat))
}


#Data are supposed to be ordered according to distance to the putative 
#pollution source
besagnewell.stat<-function(data, k)
{
	csum<-cumsum(data$Observed)-1

	#Number of regions needed to sum k cases
	l<-sum(csum<k)+1
	if(l>length(data[[1]]) ) l<-l-1

	return(c(value=l, size=l))
}
