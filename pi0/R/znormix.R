znormix=function(p, theoretical.null=TRUE, start.pi0, eps=1e-5, niter=Inf, verbose=FALSE)
{
  
    z=as.matrix(qnorm(1-p))
    z[is.infinite(z) & z<0]=min(z[is.finite(z)])
    z[is.infinite(z) & z>0]=max(z[is.finite(z)])
    G=length(z); stopifnot(G>=4)
    
    #starting values
    if(missing(start.pi0) || start.pi0<0 || start.pi0>1) {
        requireNamespace('qvalue',quitely=TRUE)
        start.pi0=qvalue(p)$pi0
    }
    if(start.pi0<=0) start.pi0=1e-3
    if(start.pi0>=1) start.pi0=1-1e-3
    if(isTRUE(theoretical.null)){
        mu1=mean(z)/(1-start.pi0)
        last.par=c(start.pi0, 0, 1, mu1, sqrt(pmax(1e-3,(var(z)-start.pi0-start.pi0*(1-start.pi0)*mu1*mu1)/(1-start.pi0))))
    }else{
        zcut=quantile(z, start.pi0)
        z0.idx=which(z<zcut)
        last.par=c(start.pi0, mean(z[z0.idx]), sd(drop(z[z0.idx])), mean(z[-z0.idx]), sd(drop(z[-z0.idx])))
    }

    #constrained EM algorithm
    iter=1
    new.par=last.par
    repeat{
#        ppee=1/(1+(1-last.par[1])/last.par[1]* 
#                exp(dnorm(z, last.par[4], last.par[5],log=TRUE)-dnorm(z, last.par[2], last.par[3],log=TRUE)))
        f0=last.par[1]*dnorm(z, last.par[2], last.par[3])
        ppee=pmin(1-1e-6, pmax(1e-6, f0/(f0+(1-last.par[1])*dnorm(z, last.par[4], last.par[5]))))
        new.par[1]=mean(ppee)
        sum.ppee=sum(ppee)
#        tmp=cov.wt(z, 1-ppee, method='ML')
#        new.par[4]=tmp$center
#        new.par[5]=sqrt(tmp$cov)
        new.par[4]=crossprod(z, 1-ppee)/(G-sum.ppee)
        new.par[5]=sqrt(crossprod((z-new.par[4])^2, 1-ppee)/(G-sum.ppee))
        if(!isTRUE(theoretical.null)){
#            tmp=cov.wt(z, ppee, method='ML')
#            new.par[2]=tmp$center
#            new.par[3]=sqrt(tmp$cov)
            new.par[2]=crossprod(z, ppee)/sum.ppee
            new.par[3]=sqrt(crossprod((z-new.par[2])^2, ppee)/sum.ppee)
            if(abs(new.par[2])>abs(new.par[4])){    # null center should be closer to 0
                tmp=new.par[2]; new.par[2]=new.par[4]; new.par[4]=tmp
                tmp=new.par[3]; new.par[3]=new.par[5]; new.par[5]=tmp
            }
        }

        if(isTRUE(verbose))
            cat('iter',iter,'\tparameters=',new.par,'\tmax.diff=',max(abs(new.par-last.par)), fill=TRUE)
        if(iter>=niter || max(abs(new.par-last.par))<eps) break
        last.par=new.par
        iter=iter+1
    }
    ord=order(ppee)
    fdr=numeric(G)
    fdr[ord]=cumsum(ppee[ord])/(1:G)
    names(new.par)=c('pi0', 'mean.z0', 'sd.z0', 'mean.z1', 'sd.z1')
    attr(new.par, 'theoretical.null')=isTRUE(theoretical.null)
    attr(new.par, 'converged')=iter<niter
    attr(new.par, 'iter')=iter
    attr(new.par, 'call')=match.call()
    attr(new.par, 'lfdr')=ppee
    attr(new.par, 'fdr')=fdr
    class(new.par)='znormix'
    new.par
}

