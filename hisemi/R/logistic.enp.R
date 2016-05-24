logistic.enp=function(log.spar, enps, maximum, minimum=2, eps=1e-8)
{
    log.spar.bak=log.spar
    enps.bak=enps
    finite.idx=is.finite(log.spar)
    log.spar=log.spar[finite.idx]
    enps=enps[finite.idx]

    n.spar=length(log.spar)
    ord=order(log.spar)
    log.spar=log.spar[ord]
    enps=pmin(pmax(enps[ord], minimum+eps), maximum-eps)

    ##########  find mode
    #require(Iso)
    uf=ufit(enps,x=log.spar)
    goodenp.idx=logical(n.spar)
    sm=which(log.spar==uf$mode)[1]; if(is.na(sm)) sm=1
    goodenp.idx[sm:n.spar]=TRUE

    n=sum(goodenp.idx)
    if(n<=1) {
#        warning("not enough effective data points. return pava results only")
        ans=enps.bak
        ans[finite.idx]=(pava(enps,decreasing=TRUE))
        attr(ans,'log.spar')=log.spar
        attr(ans,'rate')=NA
        attr(ans,'mids')=NA
        attr(ans,'pow')=NA
        attr(ans,'fit')=uf
        attr(ans,'good.idx')=goodenp.idx
        attr(ans,'mode')=uf$mode
        return(ans)
    }
    y=enps[goodenp.idx]
    x=log.spar[goodenp.idx]

    midenp=(maximum+minimum)*.5
    idx=y>midenp

    ########### find cheap start values
    fit.one.data=function(i){
        log((maximum-minimum)/(y[i]-minimum)-1)/(x[i]-mid.s)
    }
    if(any(idx) && !all(idx) ){
        big=tail(which(idx),1)
        small=big+1
    }else if (all(idx)){
        small=n
        big=n-1
    }else if (all(!idx)){
        big=1
        small=2
    }
    mid.s=(x[big]*(midenp-y[small])+x[small]*(y[big]-midenp))/(y[big]-y[small]) #linear interpolation btwn the middle two points
    rate0=median(sapply(1:n, fit.one.data),na.rm=TRUE)



    nls.fit=nls(y~minimum+(maximum-minimum)/(1+exp(rate*(x-mids)))^pow, start=list(rate=rate0, mids=mid.s, pow=1),
            lower=c(0,-Inf,0), algorithm='port')
    ans=enps.bak
    ans[finite.idx]=predict(nls.fit, newdata=list(x=log.spar))
    attr(ans,'log.spar')=log.spar
    attr(ans,'rate')=coef(nls.fit)['rate']
    attr(ans,'mids')=coef(nls.fit)['mids']
    attr(ans,'pow')=coef(nls.fit)['pow']
    attr(ans,'fit')=nls.fit
    attr(ans,'good.idx')=goodenp.idx
    attr(ans,'mode')=uf$mode
    ans
}

        