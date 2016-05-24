#
#p-value computed using a Gumbel distribution as in Jung (2009, SiM)
#

kn.gumbel.iscluster<-function(data, idx, idxorder, alpha, fractpop, use.poisson=TRUE,
        model="poisson", R, mle)
{
        #Indexes of regions in the balls ordered by distance to the centre
        localidx<-idxorder[idx[idxorder]]


        knboot<-switch(model,
        permutation = boot(data[idxorder,], statistic=kullnagar.boot, R=R, fractpop=fractpop, use.poisson=use.poisson),
        multinomial=boot(data[localidx,], statistic=kullnagar.pboot, sim="parametric", ran.gen=multinom.sim,  R=R, fractpop=fractpop, use.poisson=use.poisson, mle=list(n=mle$n, p=mle$p[localidx]) ),
        poisson = boot(data[localidx,], statistic=kullnagar.pboot, sim="parametric", ran.gen=poisson.sim,  R=R, fractpop=fractpop, use.poisson=use.poisson, mle=list(n=mle$n, lambda=mle$lambda[localidx]) ),
        negbin = boot(data[localidx,], statistic=kullnagar.pboot, sim="parametric", ran.gen=negbin.sim,  R=R, fractpop=fractpop, use.poisson=use.poisson, mle=list(n=mle$n, size=mle$size,prob=mle$prob[localidx]) )
        )

        if(is.null(knboot$t0))
                return(c(NA, NA, NA, NA))

	LLRmean<-mean(knboot$t[,1])
	LLRsd<-sd(knboot$t[,1])
	bhat<-LLRsd*sqrt(6)/pi
	ahat<-LLRmean-0.5772*bhat

        pvalue<-1-exp(-exp(-(knboot$t0[1]-ahat)/bhat))
        return(c( knboot$t0[1], alpha>pvalue, pvalue, knboot$t0[2]) )
}

