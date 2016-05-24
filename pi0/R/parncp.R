parncpt=function(tstat, df, zeromean=TRUE,  ...)
{   stopifnot(all(df>0))
    if(!zeromean && any(is.infinite(df)) && !all(is.infinite(df)) ) {
        df[is.infinite(df)]=500
    }
     method=c('L-BFGS-B')
#     method=match.arg(method)
    if       (method=='EM') {
        stop("EM algorithm not implemented")
#        parncpt.em(tstat,df,zeromean,...)
    }else if (method=='NR') {
        stop("Newton-Raphson algorithm not implemented")
#        parncpt.nr(tstat,df,zeromean,...)
    }else if (method=='L-BFGS-B') {
        if(zeromean) parncpt.bfgs.0mean(tstat,df,...) else parncpt.bfgs.non0mean(tstat,df,...)
    }
}


parncpt.bfgs.non0mean=function(tstat,df,starts, grids, approximation='int2',...)
{
    G=max(c(length(tstat),length(df)))

    dt.null=dt(tstat,df)
    obj=function(parms){
            pi0=parms[1]; mu.ncp=parms[2]; sd.ncp=parms[3]; scale.fact=sqrt(1+sd.ncp*sd.ncp)
            Lik=pi0*dt.null+(1-pi0)*dtn.mix(tstat,df,mu.ncp,sd.ncp,FALSE,approximation)
#            Lik=pi0*dt.null+(1-pi0)*dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact  
            ans=-sum(log(Lik))
            if(!is.finite(ans)){ ans=-sum(log(pi0*dt.null+(1-pi0)*dtn.mix(tstat,df,mu.ncp,sd.ncp,FALSE,approximation='none')))  }
            ans
        }
        
    deriv.non0mean=function(parms) {
        pi0=parms[1]; mu.ncp=parms[2]; sd.ncp=parms[3]; scale.fact=sqrt(1+sd.ncp*sd.ncp); s2=scale.fact*scale.fact
        dt.alt=dtn.mix(tstat, df, mu.ncp, sd.ncp, FALSE, approximation)
#        dt.alt=dt.int(tstat/scale.fact,df,mu.ncp/scale.fact)/scale.fact
        f=(pi0*dt.null+(1-pi0)*dt.alt)

        der.pi0=sum( (dt.alt-dt.null) / f )  ## correct

        if(all(is.infinite(df))){
            z.std=(tstat-mu.ncp)/scale.fact
            der.mu=-(1-pi0)*sum( dt.alt*z.std/scale.fact / f)
            der.scale=(1-pi0)*sum( dt.alt*(1-z.std*z.std)/scale.fact / f)
        }else{ df[is.infinite(df)]=500
            df.half=df/2; t2=tstat*tstat; t2vs2=t2+df*s2
            logK2=df.half*log(df.half)-.5*log(pi/2)-lgamma(df.half)
            logC=logK2-(df.half+.5)*log(t2/s2+df)- df.half*mu.ncp*mu.ncp/t2vs2
    #        integral.xv1=.C('intTruncNormVec',n=as.integer(G),r=rep(as.integer(df+1),length=G), mu=tstat*mu.ncp/scale.fact/sqrt(t2vs2),
    #                                    sd=rep(as.double(1),length=G), lower=numeric(G), upper=rep(Inf,length=G), ans=numeric(G),NAOK=TRUE)$ans
            integral.xv1=mTruncNorm.int2(r=df+1, mu=tstat*mu.ncp/scale.fact/sqrt(t2vs2),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
            
            der.mu=-sum((1-pi0)/f/s2*(tstat*exp(logC)/sqrt(t2vs2)*integral.xv1-mu.ncp*dt.alt))    ## correct

            der.scale=-sum((1-pi0)/f        /s2/scale.fact/t2vs2*(#*dhs.ds)
                dt.alt*(s2*df*(t2-s2)+mu.ncp*mu.ncp*t2vs2)-exp(logC)*mu.ncp*tstat*(t2vs2+df*s2)/sqrt(t2vs2)*integral.xv1)
            )
        }
        der.sd=sd.ncp/scale.fact*der.scale
        c(pi0=der.pi0, mu.ncp=der.mu, sd.ncp=der.sd)
    }
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, -2, 1e-3), upper=c(1-1e-3, 2, 2), ngrid=c(5,5,5))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
        starts=grid.search(obj, default.grids$lower, default.grids$upper, default.grids$ngrid)
    }
    names(starts)=c('pi0','mu.ncp','sd.ncp')
    optimFit=try(optim(starts,obj,gr=deriv.non0mean, method='L-BFGS-B',lower=c(0,-Inf,0),upper=c(1,Inf,Inf),hessian=TRUE,...))
    if(class(optimFit)=='try-error'){
        optimFit=try(nlminb(starts,obj,deriv.non0mean,lower=c(0,-Inf,0),upper=c(1,Inf,Inf), ...))
    }
    if(class(optimFit)=='try-error'){
        return(NA_real_)
    }

    ll=-optimFit$value
    attr(ll,'df')=3
	attr(ll,'nobs')=G
    class(ll)='logLik'

#library("numDeriv")
#    tmp=make.link('logit'); logit=tmp$linkfun; logitinv=tmp$linkinv; dlogitinv=tmp$mu.eta
#    obj.nobound=function(par)obj(c(logitinv(par[1]),par[2], exp(par[3])))
#    app.hess.nobound=hessian(obj.nobound, c(logit(optimFit$par[1]), optimFit$par[2], log(optimFit$par[3])))  ## need to consider hitting boundaries
#    app.hess=app.hess.nobound/tcrossprod(c(dlogitinv(logit(optimFit$par[1])), 1, optimFit$par[3]))
#
    ans=list(pi0=optimFit$par[1], mu.ncp=optimFit$par[2], sd.ncp=optimFit$par[3], data=list(tstat=tstat, df=df), 
             logLik=ll, enp=3, par=optimFit$par,
             obj=obj, gradiant=deriv.non0mean(optimFit$par), hessian=optimFit$hessian,nobs=G)
    class(ans)=c('parncpt','ncpest')
    ans
}

parncpt.bfgs.0mean=function(tstat,df, starts, grids, approximation='int2',...)
{
    G=max(c(length(tstat),length(df)))

    dt.null=dt(tstat,df)

    obj=function(parms){
                pi0=parms[1]; mu.ncp=0; sd.ncp=parms[2]; s2=1+sd.ncp*sd.ncp; scale.fact=sqrt(s2)
                Lik=pi0*dt(tstat,df)+(1-pi0)*dt(tstat/scale.fact,df)/scale.fact
                -sum(log(Lik))
    }
    deriv.0mean=function(parms){
        pi0=parms[1]; mu.ncp=0; sd.ncp=parms[2]; s2=1+sd.ncp*sd.ncp; scale.fact=sqrt(s2)
        dt.alt=dt(tstat/scale.fact,df)/scale.fact
        f=(pi0*dt.null+(1-pi0)*dt.alt)

        der.pi0=sum( (dt.alt-dt.null) / f )  ## correct

        if(all(is.infinite(df))){
            z.std=(tstat-mu.ncp)/scale.fact
#            der.mu=-(1-pi0)*sum( dt.alt*z.std/scale.fact / f)
            der.scale=(1-pi0)*sum( dt.alt*(1-z.std*z.std)/scale.fact / f)
        }else{ df[is.infinite(df)]=500
            df.half=df/2; t2=tstat*tstat; t2vs2=t2+df*s2
            logK2=df.half*log(df.half)-.5*log(pi/2)-lgamma(df.half)
            logC=logK2-(df.half+.5)*log(t2/s2+df)- df.half*mu.ncp*mu.ncp/t2vs2
    #        integral.xv1=.C('intTruncNormVec',n=as.integer(G),r=rep(as.integer(df+1),length=G), mu=tstat*mu.ncp/scale.fact/sqrt(t2vs2),
    #                                    sd=rep(as.double(1),length=G), lower=numeric(G), upper=rep(Inf,length=G), ans=numeric(G),NAOK=TRUE)$ans
            integral.xv1=mTruncNorm.int2(r=df+1, mu=tstat*mu.ncp/scale.fact/sqrt(t2vs2),
                                        sd=1, lower=0, upper=Inf, takeLog=TRUE, ndiv=8)
    #        der.mu=-sum((1-pi0)/f/s2*(tstat*exp(logC)/sqrt(t2vs2)*integral.xv1-mu.ncp*dt.alt))    ## correct

            der.scale=-sum((1-pi0)/f        /s2/scale.fact/t2vs2*(#*dhs.ds)
                dt.alt*(s2*df*(t2-s2)+mu.ncp*mu.ncp*t2vs2)-exp(logC)*mu.ncp*tstat*(t2vs2+df*s2)/sqrt(t2vs2)*integral.xv1)
            )
        }
        der.sd=sd.ncp/scale.fact*der.scale
        c(der.pi0, der.sd)
    }
    if(missing(starts)) {
        default.grids=list(lower=c(1e-3, 1e-3), upper=c(1-1e-3, 2), ngrid=c(15,15))
        if(!missing(grids)) for(nn in names(grids)) default.grids[[nn]]=grids[[nn]]
        starts=grid.search(obj, default.grids$lower, default.grids$upper, default.grids$ngrid)
    }
    names(starts)=c('pi0','sd.ncp')
    optimFit=optim(starts,obj,gr=deriv.0mean, method='L-BFGS-B', lower=c(0,0),upper=c(1,Inf), hessian=TRUE,...)

    ll=-optimFit$value
    attr(ll,'df')=2
	attr(ll,'nobs')=G
    class(ll)='logLik'
    
#library("numDeriv")
#    tmp=make.link('logit'); logit=tmp$linkfun; logitinv=tmp$linkinv; dlogitinv=tmp$mu.eta
#    obj.nobound=function(par)obj(c(logitinv(par[1]),exp(par[2])))
#    app.hess.nobound=hessian(obj.nobound, c(logit(optimFit$par[1]), log(optimFit$par[2])))  ## need to consider hitting boundaries
#    app.hess=app.hess.nobound/tcrossprod(c(dlogitinv(logit(optimFit$par[1])), optimFit$par[2]))
#
    ans=list(pi0=optimFit$par[1], mu.ncp=0, sd.ncp=optimFit$par[2], data=list(tstat=tstat, df=df), 
             logLik=ll, enp=2, par=optimFit$par,
             obj=obj, gradiant=deriv.0mean(optimFit$par), hessian=optimFit$hessian,nobs=G)
    class(ans)=c('parncpt','ncpest')
    ans
}

## vcov and logLik are defined for the ncpest class in nparncp.R
#vcov.parncp=function(obj)
#{
#    obj$hessian
#}
#logLik.parncp=function(obj)
#{
#    obj$logLik
#}
#coef.parncp=#coefficients.parncp=
#function(obj)
#{
#    obj$par
#}
fitted.parncpt=#fitted.values.parncpt=
function(object, ...)
{
    object$pi0*dt(object$data$tstat, object$data$df)+(1-object$pi0)*dtn.mix(object$data$tstat, object$data$df,object$mu.ncp,object$sd.ncp,FALSE,...)
}
#lfdr=ppee=function(object)UseMethod('lfdr')
#lfdr.parncpt=ppee.parncpt=function(object)
#{
#    pmin(pmax(object$pi0*dt(object$data$tstat, object$data$df)/fitted(object), 0), 1)
#}
summary.parncpt=function(object,...)
{
    cat("pi0 (proportion of null hypotheses) =", object$pi0, fill=TRUE)
    cat("mu.ncp (mean of noncentrality parameters) =", object$mu.ncp, fill=TRUE)
    cat("sd.ncp (SD of noncentrality parameters) =", object$sd.ncp, fill=TRUE)
    invisible(object)
}
print.parncpt=function(x,...)
{
    summary.parncpt(x,...)
}
plot.parncpt=function(x,...)
{
#    dev.new(width=8, height=4)
    op=par(mfrow=c(1,2))
    hist(x$data$tstat, pr=TRUE, br=min(c(max(c(20, length(x$data$tstat)/100)), 200)), xlab='t',main='t-statistics')
    ord=order(x$data$tstat)
    lines(x$data$tstat[ord], fitted.parncpt(x)[ord], col='red', lwd=2)
    d.ncp=function(d) dnorm(d, x$mu.ncp, x$sd.ncp)
    curve(d.ncp, min(x$data$tstat), max(x$data$tstat), 500, xlab=expression(delta), ylab='density',main='noncentrality parameters')
    abline(v=c(0, x$mu.ncp), lty=1:2)
    par(op)
    invisible(x)
}

grid.search=function(obj, lower, upper, ngrid, ...)
{
    p=max(c(length(lower),length(upper),length(ngrid)))
    lower=rep(lower, length=p)
    upper=rep(upper, length=p)
    ngrid=rep(ngrid, length=p)

    knot.list=list()
    for(i in 1:p) knot.list[[i]]=seq(from=lower[i], to=upper[i], length=ngrid[i])
    names(knot.list)=paste(names(lower),names(upper),names(ngrid),sep='.')
    
    grids=do.call(expand.grid, knot.list)
    ans=apply(grids, 1, obj, ...)
    if(sum(ans==min(ans))>1) warning('multiple minimums found in grid search')
    return(unlist(grids[which.min(ans),]))
}

