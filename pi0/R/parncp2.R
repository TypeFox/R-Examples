parncpt2=function(tstat, df, common=c('mean','sd'), ...)
{   stopifnot(all(df>0))
    # if(any(is.infinite(df)) && !all(is.infinite(df)) ) {
        # df[is.infinite(df)]=500
    # }
     method=c('constrOptim')
	 if(length(common)>0L){
		common=tolower(common)
		common=match.arg(common, c('mean','sd','none'), several.ok=TRUE)
		if('none'%in%common && length(common)>1L) stop('"common" argument incompatible.')
		common=sort(common)
		if('none'%in%common)common=NULL
	 }
#     method=match.arg(method)
    if       (method=='EM') {
        stop("EM algorithm not implemented")
#        parncpt.em(tstat,df,zeromean,...)
    }else if (method=='NR') {
        stop("Newton-Raphson algorithm not implemented")
#        parncpt.nr(tstat,df,zeromean,...)
    }else if (method=='constrOptim') {
		if(is.null(common)){
			parncpt2.constrOptim(tstat,df,...)
		}else if(all(common==c('mean','sd'))) {
			parncpt2.constrOptim.CMD(tstat,df,...) 
		}else if(common=='mean'){
			parncpt2.constrOptim.CM(tstat,df,...) 
		}else if(common=='sd'){
			parncpt2.constrOptim.CD(tstat,df,...) 
		}else stop('something wrong with "common" arg.')
    }
}

parncpt2.constrOptim.CM=function(tstat,df,starts, grids, approximation='int2',...)
{
    G=max(c(length(tstat),length(df)))

    dt.null=dt(tstat,df)
    obj=function(parms){
            pi0=parms[1]; pi1=parms[2]; mu1.ncp=parms[3]; sd1.ncp=parms[4]; mu2.ncp=-parms[3]; sd2.ncp=parms[5]; 
			 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
            Lik=pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation)
#            Lik=pi0*dt.null+(1-pi0)*dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation='none')+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.obj=function(parms) {names(parms)=NULL
       pi0=parms[[1]]; pi1=parms[[2]]; mu1.ncp=parms[[3]]; sd1.ncp=parms[[4]]; mu2.ncp=-parms[[3]]; sd2.ncp=parms[[5]]; 
		 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
        s2_1=scale.fact1*scale.fact1; s2_2=scale.fact2*scale.fact2
        dt.alt1=dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)
        dt.alt2=dtn.mix(tstat, df, mu2.ncp, sd2.ncp, FALSE, approximation)
#        dt.alt=dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact
        f=(pi0*dt.null+pi1*dt.alt1+(1-pi0-pi1)*dt.alt2)

        der.pi0=sum( (dt.alt2-dt.null) / f )  ## correct
        der.pi1=sum( (dt.alt2-dt.alt1) / f )  ## correct

        if(all(is.infinite(df))){
            z.std1=(tstat-mu1.ncp)/scale.fact1
            der.mu1=-pi1*sum( dt.alt1*z.std1/scale.fact1 / f)
            der.scale1=pi1*sum( dt.alt1*(1-z.std1*z.std1)/scale.fact1 / f)
			
            z.std2=(tstat-mu2.ncp)/scale.fact2
            der.mu2=-(1-pi0-pi1)*sum( dt.alt2*z.std2/scale.fact2 / f)
            der.scale2=(1-pi0-pi1)*sum( dt.alt2*(1-z.std2*z.std2)/scale.fact2 / f)
        }else{ df[is.infinite(df)]=500
            df.half=df/2; t2=tstat*tstat; 
            logK2=df.half*log(df.half)-.5*log(pi/2)-lgamma(df.half)
			
			 t2vs2_1=t2+df*s2_1
            logC1=logK2-(df.half+.5)*log(t2/s2_1+df)- df.half*mu1.ncp*mu1.ncp/t2vs2_1
            integral.xv1=mTruncNorm.int2(r=df+1, mu=tstat*mu1.ncp/scale.fact1/sqrt(t2vs2_1),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu1=-sum(pi1/f/s2_1*(tstat*exp(logC1)/sqrt(t2vs2_1)*integral.xv1-mu1.ncp*dt.alt1))    ## correct
            der.scale1=-sum(pi1/f        /s2_1/scale.fact1/t2vs2_1*(#*dhs.ds)
                dt.alt1*(s2_1*df*(t2-s2_1)+mu1.ncp*mu1.ncp*t2vs2_1)-exp(logC1)*mu1.ncp*tstat*(t2vs2_1+df*s2_1)/sqrt(t2vs2_1)*integral.xv1)
            )
			
			
			 t2vs2_2=t2+df*s2_2
            logC2=logK2-(df.half+.5)*log(t2/s2_2+df)- df.half*mu2.ncp*mu2.ncp/t2vs2_2
            integral.xv2=mTruncNorm.int2(r=df+1, mu=tstat*mu2.ncp/scale.fact2/sqrt(t2vs2_2),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu2=-sum((1-pi0-pi1)/f/s2_2*(tstat*exp(logC2)/sqrt(t2vs2_2)*integral.xv2-mu2.ncp*dt.alt2))    ## correct
            der.scale2=-sum((1-pi0-pi1)/f        /s2_2/scale.fact2/t2vs2_2*(#*dhs.ds)
                dt.alt2*(s2_2*df*(t2-s2_2)+mu2.ncp*mu2.ncp*t2vs2_2)-exp(logC2)*mu2.ncp*tstat*(t2vs2_2+df*s2_2)/sqrt(t2vs2_2)*integral.xv2)
            )
		}
        der.sd1=sd1.ncp/scale.fact1*der.scale1
        der.sd2=sd2.ncp/scale.fact2*der.scale2
		
        c(pi0=der.pi0, pi1=der.pi1, mu1.ncp=der.mu1-der.mu2, sd1.ncp=der.sd1, sd2.ncp=der.sd2)
    }
	#deriv.obj=function(parms)numDeriv::grad(obj,parms)
	
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3, -3, 1e-3), upper=c(1-1e-3, 1-1e-3, -1e-3, 2), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
		obj.restricted=function(parms){
			if(sum(parms[1:2])>=1) Inf else
				obj(c(parms[1],parms[2],parms[3],parms[4],parms[3]*-1,parms[4]))
		}
        starts=suppressWarnings(grid.search(obj.restricted, default.grids$lower, default.grids$upper, default.grids$ngrid))
		starts=c(starts[1:4], starts[4])
    }
	ui=rbind(diag(c(1,1,-1,1,1)), rep(-1:0, 2:3))
	ci=rep(0,6); ci[6]=-1
    names(starts)=c('pi0','pi1','mu1.ncp','sd1.ncp','sd2.ncp')
    optimFit=try(constrOptim(starts,obj,grad=deriv.obj, ui=ui,ci=ci,hessian=FALSE,...))
    # if(class(optimFit)=='try-error'){
        # optimFit=try(nlminb(starts,obj,deriv.non0mean,lower=c(0,-Inf,0),upper=c(1,Inf,Inf), ...))
    # }
    if(class(optimFit)=='try-error'){
        return(NA_real_)
    }
	optimFit$hessian=numDeriv::hessian(obj, optimFit$par)

    ll=-optimFit$value
    attr(ll,'df')=5
	attr(ll,'nobs')=G
    class(ll)='logLik'

	tau=optimFit$par[2L]/(1-optimFit$par[1L])
	mu=tau*optimFit$par[3L]-(1-tau)*optimFit$par[3L]
	Var=tau*((optimFit$par[3L]-mu)^2+optimFit$par[4L]^2)+
	(1-tau)*((-optimFit$par[3L]-mu)^2+optimFit$par[5L]^2)
	
    ans=list(pi0=optimFit$par[1], mu.ncp=mu, sd.ncp=sqrt(Var), tau.ncp=tau, pi1=optimFit$par[2], mu1.ncp=optimFit$par[3], sd1.ncp=optimFit$par[4], 
		mu2.ncp=-optimFit$par[3], sd2.ncp=optimFit$par[5], 
		data=list(tstat=tstat, df=df), 
             logLik=ll, enp=5, par=optimFit$par,
             obj=obj, gradiant=deriv.obj(optimFit$par), hessian=optimFit$hessian,nobs=G)
    class(ans)=c('parncpt2','parncpt','ncpest')
    ans
}

parncpt2.constrOptim.CD=function(tstat,df,starts, grids, approximation='int2',...)
{
    G=max(c(length(tstat),length(df)))

    dt.null=dt(tstat,df)
    obj=function(parms){
            pi0=parms[1]; pi1=parms[2]; mu1.ncp=parms[3]; sd1.ncp=parms[4]; mu2.ncp=parms[5]; sd2.ncp=parms[4]; 
			 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
            Lik=pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation)
#            Lik=pi0*dt.null+(1-pi0)*dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation='none')+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.obj=function(parms) {
        pi0=parms[[1]]; pi1=parms[[2]]; mu1.ncp=parms[[3]]; sd1.ncp=parms[[4]]; mu2.ncp=parms[[5]]; sd2.ncp=parms[[4]]; 
		 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
        s2_1=scale.fact1*scale.fact1; s2_2=scale.fact2*scale.fact2
        dt.alt1=dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)
        dt.alt2=dtn.mix(tstat, df, mu2.ncp, sd2.ncp, FALSE, approximation)
#        dt.alt=dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact
        f=(pi0*dt.null+pi1*dt.alt1+(1-pi0-pi1)*dt.alt2)

        der.pi0=sum( (dt.alt2-dt.null) / f )  ## correct
        der.pi1=sum( (dt.alt2-dt.alt1) / f )  ## correct

        if(all(is.infinite(df))){
            z.std1=(tstat-mu1.ncp)/scale.fact1
            der.mu1=-pi1*sum( dt.alt1*z.std1/scale.fact1 / f)
            der.scale1=pi1*sum( dt.alt1*(1-z.std1*z.std1)/scale.fact1 / f)
			
            z.std2=(tstat-mu2.ncp)/scale.fact2
            der.mu2=-(1-pi0-pi1)*sum( dt.alt2*z.std2/scale.fact2 / f)
            der.scale2=(1-pi0-pi1)*sum( dt.alt2*(1-z.std2*z.std2)/scale.fact2 / f)
        }else{ df[is.infinite(df)]=500
            df.half=df/2; t2=tstat*tstat; 
            logK2=df.half*log(df.half)-.5*log(pi/2)-lgamma(df.half)
			
			 t2vs2_1=t2+df*s2_1
            logC1=logK2-(df.half+.5)*log(t2/s2_1+df)- df.half*mu1.ncp*mu1.ncp/t2vs2_1
            integral.xv1=mTruncNorm.int2(r=df+1, mu=tstat*mu1.ncp/scale.fact1/sqrt(t2vs2_1),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu1=-sum(pi1/f/s2_1*(tstat*exp(logC1)/sqrt(t2vs2_1)*integral.xv1-mu1.ncp*dt.alt1))    ## correct
            der.scale1=-sum(pi1/f        /s2_1/scale.fact1/t2vs2_1*(#*dhs.ds)
                dt.alt1*(s2_1*df*(t2-s2_1)+mu1.ncp*mu1.ncp*t2vs2_1)-exp(logC1)*mu1.ncp*tstat*(t2vs2_1+df*s2_1)/sqrt(t2vs2_1)*integral.xv1)
            )
			
			
			 t2vs2_2=t2+df*s2_2
            logC2=logK2-(df.half+.5)*log(t2/s2_2+df)- df.half*mu2.ncp*mu2.ncp/t2vs2_2
            integral.xv2=mTruncNorm.int2(r=df+1, mu=tstat*mu2.ncp/scale.fact2/sqrt(t2vs2_2),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu2=-sum((1-pi0-pi1)/f/s2_2*(tstat*exp(logC2)/sqrt(t2vs2_2)*integral.xv2-mu2.ncp*dt.alt2))    ## correct
            der.scale2=-sum((1-pi0-pi1)/f        /s2_2/scale.fact2/t2vs2_2*(#*dhs.ds)
                dt.alt2*(s2_2*df*(t2-s2_2)+mu2.ncp*mu2.ncp*t2vs2_2)-exp(logC2)*mu2.ncp*tstat*(t2vs2_2+df*s2_2)/sqrt(t2vs2_2)*integral.xv2)
            )
		}
        der.sd1=sd1.ncp/scale.fact1*der.scale1
        der.sd2=sd2.ncp/scale.fact2*der.scale2
		
        c(pi0=der.pi0, pi1=der.pi1, mu1.ncp=der.mu1, sd1.ncp=der.sd1+der.sd2, mu2.ncp=der.mu2)
    }
	#deriv.obj=function(parms)numDeriv::grad(obj, parms)
	
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3, -2, 1e-3), upper=c(1-1e-3, 1-1e-3, 2, 2), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
		obj.restricted=function(parms){
			if(sum(parms[1:2])>=1) Inf else
				obj(c(parms[1],parms[2],parms[3],parms[4],parms[3]*-1,parms[4]))
		}
        starts=suppressWarnings(grid.search(obj.restricted, default.grids$lower, default.grids$upper, default.grids$ngrid))
		if(starts[3]>0){
			starts[2]=1-starts[1]-starts[2]; starts[3]=-starts[3]
		}else if(starts[3]==0) starts[3]=-1e-3
		starts=c(starts[1:4], starts[3]*-1)
    }
	ui=rbind(diag(1,5)[c(1,2,4),], rep(-1:0, 2:3), c(0,0,-1,0,1))
	ci=rep(0,5); ci[4]=-1
    names(starts)=c('pi0','pi1','mu1.ncp','sd1.ncp','mu2.ncp')
    optimFit=try(constrOptim(starts,obj,grad=deriv.obj, ui=ui,ci=ci,hessian=FALSE,...))
    # if(class(optimFit)=='try-error'){
        # optimFit=try(nlminb(starts,obj,deriv.non0mean,lower=c(0,-Inf,0),upper=c(1,Inf,Inf), ...))
    # }
    if(class(optimFit)=='try-error'){
        return(NA_real_)
    }
	optimFit$hessian=numDeriv::hessian(obj, optimFit$par)

    ll=-optimFit$value
    attr(ll,'df')=5
	attr(ll,'nobs')=G
    class(ll)='logLik'

	tau=optimFit$par[2L]/(1-optimFit$par[1L])
	mu=tau*optimFit$par[3L]+(1-tau)*optimFit$par[5L]
	Var=tau*((optimFit$par[3L]-mu)^2+optimFit$par[4L]^2)+
	(1-tau)*((optimFit$par[5L]-mu)^2+optimFit$par[4L]^2)
	
    ans=list(pi0=optimFit$par[1], mu.ncp=mu, sd.ncp=sqrt(Var), tau.ncp=tau, pi1=optimFit$par[2], mu1.ncp=optimFit$par[3], sd1.ncp=optimFit$par[4], 
		mu2.ncp=optimFit$par[5], sd2.ncp=optimFit$par[4], 
		data=list(tstat=tstat, df=df), 
             logLik=ll, enp=5, par=optimFit$par,
             obj=obj, gradiant=deriv.obj(optimFit$par), hessian=optimFit$hessian,nobs=G)
    class(ans)=c('parncpt2','parncpt','ncpest')
    ans
}



parncpt2.constrOptim.CMD=function(tstat,df,starts, grids, approximation='int2',...)
{
    G=max(c(length(tstat),length(df)))

    dt.null=dt(tstat,df)
    obj=function(parms){
            pi0=parms[1]; pi1=parms[2]; mu1.ncp=parms[3]; sd1.ncp=parms[4]; mu2.ncp=-parms[3]; sd2.ncp=parms[4]; 
			 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
            Lik=pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation)
#            Lik=pi0*dt.null+(1-pi0)*dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation='none')+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.obj=function(parms) {
        pi0=parms[[1]]; pi1=parms[[2]]; mu1.ncp=parms[[3]]; sd1.ncp=parms[[4]]; mu2.ncp=-parms[[3]]; sd2.ncp=parms[[4]]; 
		 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
        s2_1=scale.fact1*scale.fact1; s2_2=scale.fact2*scale.fact2
        dt.alt1=dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)
        dt.alt2=dtn.mix(tstat, df, mu2.ncp, sd2.ncp, FALSE, approximation)
#        dt.alt=dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact
        f=(pi0*dt.null+pi1*dt.alt1+(1-pi0-pi1)*dt.alt2)

        der.pi0=sum( (dt.alt2-dt.null) / f )  ## correct
        der.pi1=sum( (dt.alt2-dt.alt1) / f )  ## correct

        if(all(is.infinite(df))){
            z.std1=(tstat-mu1.ncp)/scale.fact1
            der.mu1=-pi1*sum( dt.alt1*z.std1/scale.fact1 / f)
            der.scale1=pi1*sum( dt.alt1*(1-z.std1*z.std1)/scale.fact1 / f)
			
            z.std2=(tstat-mu2.ncp)/scale.fact2
            der.mu2=-(1-pi0-pi1)*sum( dt.alt2*z.std2/scale.fact2 / f)
            der.scale2=(1-pi0-pi1)*sum( dt.alt2*(1-z.std2*z.std2)/scale.fact2 / f)
        }else{ df[is.infinite(df)]=500
            df.half=df/2; t2=tstat*tstat; 
            logK2=df.half*log(df.half)-.5*log(pi/2)-lgamma(df.half)
			
			 t2vs2_1=t2+df*s2_1
            logC1=logK2-(df.half+.5)*log(t2/s2_1+df)- df.half*mu1.ncp*mu1.ncp/t2vs2_1
            integral.xv1=mTruncNorm.int2(r=df+1, mu=tstat*mu1.ncp/scale.fact1/sqrt(t2vs2_1),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu1=-sum(pi1/f/s2_1*(tstat*exp(logC1)/sqrt(t2vs2_1)*integral.xv1-mu1.ncp*dt.alt1))    ## correct
            der.scale1=-sum(pi1/f        /s2_1/scale.fact1/t2vs2_1*(#*dhs.ds)
                dt.alt1*(s2_1*df*(t2-s2_1)+mu1.ncp*mu1.ncp*t2vs2_1)-exp(logC1)*mu1.ncp*tstat*(t2vs2_1+df*s2_1)/sqrt(t2vs2_1)*integral.xv1)
            )
			
			
			 t2vs2_2=t2+df*s2_2
            logC2=logK2-(df.half+.5)*log(t2/s2_2+df)- df.half*mu2.ncp*mu2.ncp/t2vs2_2
            integral.xv2=mTruncNorm.int2(r=df+1, mu=tstat*mu2.ncp/scale.fact2/sqrt(t2vs2_2),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu2=-sum((1-pi0-pi1)/f/s2_2*(tstat*exp(logC2)/sqrt(t2vs2_2)*integral.xv2-mu2.ncp*dt.alt2))    ## correct
            der.scale2=-sum((1-pi0-pi1)/f        /s2_2/scale.fact2/t2vs2_2*(#*dhs.ds)
                dt.alt2*(s2_2*df*(t2-s2_2)+mu2.ncp*mu2.ncp*t2vs2_2)-exp(logC2)*mu2.ncp*tstat*(t2vs2_2+df*s2_2)/sqrt(t2vs2_2)*integral.xv2)
            )
		}
        der.sd1=sd1.ncp/scale.fact1*der.scale1
        der.sd2=sd2.ncp/scale.fact2*der.scale2
		
        c(pi0=der.pi0, pi1=der.pi1, mu1.ncp=der.mu1-der.mu2, sd1.ncp=der.sd1+der.sd2)
    }
	#deriv.obj=function(parms)numDeriv::grad(obj, parms)
	
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3, -3, 1e-3), upper=c(1-1e-3, 1-1e-3, -1e-3, 2), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
		 tmp=function(parms)if(sum(parms[1:2])>=1) Inf else obj(parms)
        starts=suppressWarnings(grid.search(tmp, default.grids$lower, default.grids$upper, default.grids$ngrid))
    }
	ui=rbind(diag(1,4), rep(-1:0, each=2)); ui[3L,3L]=-1
	ci=rep(0,5); ci[5]=-1
    names(starts)=c('pi0','pi1','mu1.ncp','sd1.ncp')
    optimFit=try(constrOptim(starts,obj,grad=deriv.obj, ui=ui,ci=ci,hessian=FALSE,...))
    # if(class(optimFit)=='try-error'){
        # optimFit=try(nlminb(starts,obj,deriv.non0mean,lower=c(0,-Inf,0),upper=c(1,Inf,Inf), ...))
    # }
    if(class(optimFit)=='try-error'){
        return(NA_real_)
    }
	optimFit$hessian=numDeriv::hessian(obj, optimFit$par)

    ll=-optimFit$value
    attr(ll,'df')=4
	attr(ll,'nobs')=G
    class(ll)='logLik'

	tau=optimFit$par[2L]/(1-optimFit$par[1L])
	mu=tau*optimFit$par[3L]-(1-tau)*optimFit$par[3L]
	Var=tau*((optimFit$par[3L]-mu)^2+optimFit$par[4L]^2)+
	(1-tau)*((-optimFit$par[3L]-mu)^2+optimFit$par[4L]^2)
	
    ans=list(pi0=optimFit$par[1], mu.ncp=mu, sd.ncp=sqrt(Var), tau.ncp=tau, pi1=optimFit$par[2], mu1.ncp=optimFit$par[3], sd1.ncp=optimFit$par[4], 
		mu2.ncp=-optimFit$par[3], sd2.ncp=optimFit$par[4], 
		data=list(tstat=tstat, df=df), 
             logLik=ll, enp=4, par=optimFit$par,
             obj=obj, gradiant=deriv.obj(optimFit$par), hessian=optimFit$hessian,nobs=G)
    class(ans)=c('parncpt2','parncpt','ncpest')
    ans
}

parncpt2.constrOptim=function(tstat,df,starts, grids, approximation='int2',...)
{
    G=max(c(length(tstat),length(df)))

    dt.null=dt(tstat,df)
    obj=function(parms){
            pi0=parms[1]; pi1=parms[2]; mu1.ncp=parms[3]; sd1.ncp=parms[4]; mu2.ncp=parms[5]; sd2.ncp=parms[6]; 
			 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
            Lik=pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation)
#            Lik=pi0*dt.null+(1-pi0)*dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dt.null+pi1*dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation='none')+(1-pi0-pi1)*dtn.mix(tstat,df,mu2.ncp,sd2.ncp,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.obj=function(parms) {
		pi0=parms[[1]]; pi1=parms[[2]]; mu1.ncp=parms[[3]]; sd1.ncp=parms[[4]]; mu2.ncp=parms[[5]]; sd2.ncp=parms[[6]]; 
		 scale.fact1=sqrt(1+sd1.ncp*sd1.ncp); scale.fact2=sqrt(1+sd2.ncp*sd2.ncp); 
        s2_1=scale.fact1*scale.fact1; s2_2=scale.fact2*scale.fact2
        dt.alt1=dtn.mix(tstat, df, mu1.ncp, sd1.ncp, FALSE, approximation)
        dt.alt2=dtn.mix(tstat, df, mu2.ncp, sd2.ncp, FALSE, approximation)
#        dt.alt=dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact
        f=(pi0*dt.null+pi1*dt.alt1+(1-pi0-pi1)*dt.alt2)

        der.pi0=sum( (dt.alt2-dt.null) / f )  ## correct
        der.pi1=sum( (dt.alt2-dt.alt1) / f )  ## correct

        if(all(is.infinite(df))){
            z.std1=(tstat-mu1.ncp)/scale.fact1
            der.mu1=-pi1*sum( dt.alt1*z.std1/scale.fact1 / f)
            der.scale1=pi1*sum( dt.alt1*(1-z.std1*z.std1)/scale.fact1 / f)
			
            z.std2=(tstat-mu2.ncp)/scale.fact2
            der.mu2=-(1-pi0-pi1)*sum( dt.alt2*z.std2/scale.fact2 / f)
            der.scale2=(1-pi0-pi1)*sum( dt.alt2*(1-z.std2*z.std2)/scale.fact2 / f)
        }else{ df[is.infinite(df)]=500
            df.half=df/2; t2=tstat*tstat; 
            logK2=df.half*log(df.half)-.5*log(pi/2)-lgamma(df.half)
			
			 t2vs2_1=t2+df*s2_1
            logC1=logK2-(df.half+.5)*log(t2/s2_1+df)- df.half*mu1.ncp*mu1.ncp/t2vs2_1
            integral.xv1=mTruncNorm.int2(r=df+1, mu=tstat*mu1.ncp/scale.fact1/sqrt(t2vs2_1),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu1=-sum(pi1/f/s2_1*(tstat*exp(logC1)/sqrt(t2vs2_1)*integral.xv1-mu1.ncp*dt.alt1))    ## correct
            der.scale1=-sum(pi1/f        /s2_1/scale.fact1/t2vs2_1*(#*dhs.ds)
                dt.alt1*(s2_1*df*(t2-s2_1)+mu1.ncp*mu1.ncp*t2vs2_1)-exp(logC1)*mu1.ncp*tstat*(t2vs2_1+df*s2_1)/sqrt(t2vs2_1)*integral.xv1)
            )
			
			
			 t2vs2_2=t2+df*s2_2
            logC2=logK2-(df.half+.5)*log(t2/s2_2+df)- df.half*mu2.ncp*mu2.ncp/t2vs2_2
            integral.xv2=mTruncNorm.int2(r=df+1, mu=tstat*mu2.ncp/scale.fact2/sqrt(t2vs2_2),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            der.mu2=-sum((1-pi0-pi1)/f/s2_2*(tstat*exp(logC2)/sqrt(t2vs2_2)*integral.xv2-mu2.ncp*dt.alt2))    ## correct
            der.scale2=-sum((1-pi0-pi1)/f        /s2_2/scale.fact2/t2vs2_2*(#*dhs.ds)
                dt.alt2*(s2_2*df*(t2-s2_2)+mu2.ncp*mu2.ncp*t2vs2_2)-exp(logC2)*mu2.ncp*tstat*(t2vs2_2+df*s2_2)/sqrt(t2vs2_2)*integral.xv2)
            )
		}
        der.sd1=sd1.ncp/scale.fact1*der.scale1
        der.sd2=sd2.ncp/scale.fact2*der.scale2
		
        c(pi0=der.pi0, pi1=der.pi1, mu1.ncp=der.mu1, sd1.ncp=der.sd1, mu2.ncp=der.mu2, sd2.ncp=der.sd2)
    }
	
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3, -2, 1e-3), upper=c(1-1e-3, 1-1e-3, 2, 2), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
		obj.restricted=function(parms){
			if(sum(parms[1:2])>1) Inf else
				obj(c(parms[1],parms[2],parms[3],parms[4],parms[3]*-1,parms[4]))
		}
        starts=suppressWarnings(grid.search(obj.restricted, default.grids$lower, default.grids$upper, default.grids$ngrid))
		if(sum(starts[1:2])==1) starts[1:2]=starts[1:2]*.999
		if(starts[3]>0){
			starts[2]=1-starts[1]-starts[2]; starts[3]=-starts[3]
		}else if(starts[3]==0) starts[3]=-1e-3
		starts=c(starts[1:4], starts[3]*-1, starts[4])
    }
	ui=rbind(diag(1,6)[-c(3,5),], rep(-1:0, c(2,4)), c(0,0,-1,0,1,0))
	ci=rep(0,6); ci[5]=-1
    names(starts)=c('pi0','pi1','mu1.ncp','sd1.ncp','mu2.ncp','sd2.ncp')
    optimFit=try(constrOptim(starts,obj,grad=deriv.obj, ui=ui,ci=ci,hessian=FALSE,...))
    # if(class(optimFit)=='try-error'){
        # optimFit=try(nlminb(starts,obj,deriv.non0mean,lower=c(0,-Inf,0),upper=c(1,Inf,Inf), ...))
    # }
    if(class(optimFit)=='try-error'){
        return(NA_real_)
    }
	optimFit$hessian=numDeriv::hessian(obj, optimFit$par)

    ll=-optimFit$value
    attr(ll,'df')=6
	attr(ll,'nobs')=G
    class(ll)='logLik'

	tau=optimFit$par[2L]/(1-optimFit$par[1L])
	mu=tau*optimFit$par[3L]+(1-tau)*optimFit$par[5L]
	Var=tau*((optimFit$par[3L]-mu)^2+optimFit$par[4L]^2)+
	(1-tau)*((optimFit$par[5L]-mu)^2+optimFit$par[6L]^2)
	
    ans=list(pi0=optimFit$par[1], mu.ncp=mu, sd.ncp=sqrt(Var), tau.ncp=tau, pi1=optimFit$par[2], mu1.ncp=optimFit$par[3], sd1.ncp=optimFit$par[4], 
		mu2.ncp=optimFit$par[5], sd2.ncp=optimFit$par[6], 
		data=list(tstat=tstat, df=df), 
             logLik=ll, enp=6, par=optimFit$par,
             obj=obj, gradiant=deriv.obj(optimFit$par), hessian=optimFit$hessian,nobs=G)
    class(ans)=c('parncpt2','parncpt','ncpest')
    ans
}

fitted.parncpt2=#fitted.values.parncpt=
function(object, ...)
{
    object$pi0*dt(object$data$tstat, object$data$df)+
    object$pi1*dtn.mix(object$data$tstat, object$data$df,object$mu1.ncp,object$sd1.ncp,FALSE,...)+
	(1-object$pi0-object$pi1)*dtn.mix(object$data$tstat, object$data$df,object$mu2.ncp,object$sd2.ncp,FALSE,...)
}

summary.parncpt2=function(object,...)
{
    cat("pi0 (proportion of null hypotheses) =", object$pi0, "\n")
    cat("mu.ncp (mean of noncentrality parameters) =", object$mu.ncp, "\n")
    cat("sd.ncp (SD of noncentrality parameters) =", object$sd.ncp, "\n")
	 cat(sprintf("tau.ncp=%.5f; mu1.ncp=%.3f; sd1.ncp=%.2f; mu2.ncp=%.3f; sd2.ncp=%.2f",object$tau,object$mu1.ncp,object$sd1.ncp,object$mu2.ncp,object$sd2.ncp),"\n")
    invisible(object)
}
print.parncpt2=function(x,...)
{
    summary.parncpt2(x,...)
}
plot.parncpt2=function(x,...)
{
	tau=x$tau.ncp
	mu=x$mu.ncp
#    dev.new(width=8, height=4)
    op=par(mfrow=c(1,2))
    hist(x$data$tstat, pr=TRUE, br=min(c(max(c(20, length(x$data$tstat)/100)), 200)), xlab='t',main='t-statistics')
    ord=order(x$data$tstat)
    lines(x$data$tstat[ord], fitted.parncpt2(x)[ord], col='red', lwd=2)
    d.ncp=function(d) tau*dnorm(d, x$mu1.ncp, x$sd1.ncp)+(1-tau)*dnorm(d, x$mu2.ncp, x$sd2.ncp)
    curve(d.ncp, min(x$data$tstat), max(x$data$tstat), 500, xlab=expression(delta), ylab='density',main='noncentrality parameters')
    abline(v=c(0, mu), lty=1:2)
    par(op)
    invisible(x)
}
