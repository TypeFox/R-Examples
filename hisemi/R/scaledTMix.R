
scaledTMix.null=function(tstat,df
    ,starts=list(pi0=seq(.1,.99,length=20), scale=2^seq(.01, log2(max(abs(tstat))),length=20))
    )
{

    #### initializing central t and non-central t
    dt0=dt(tstat,df)
    t2=tstat*tstat
    G=length(t2)

    logLik=function(parms)
    {
        pi0=parms[1]; 
        if(length(parms)>1){
            scale.fact=parms[2]; #rest.parms=if(length(parms)>1)parms[-1] else NULL
            -2*sum(log(
                pi0*dt0+(1-pi0)*dt(tstat/scale.fact,df)/scale.fact
            #                rest1.delta%*%(pncp(deltas[-1],rest.parms)-pncp(deltas[-length(deltas)],rest.parms))
            ))
        }else{      ## this relies on pre-computed dt1 matrix
            (-2*(colSums(log(pi0/(1-pi0)*dt0+dt1))+G*log(1-pi0)))
        }
    }



    if(is.null(names(starts))) names(starts)=c('pi0','scale')

    #cat("searching for good starts...")        
        scale.grids=starts[['scale']]
        tstat.scale=outer(tstat,scale.grids,'/')
        dt1=sweep(dt(tstat.scale,df),2,scale.grids,'/')
        
        pi0.grids=starts[['pi0']]
        starts.rslt=sapply(pi0.grids,logLik)
        dim(starts.rslt)=c(length(scale.grids),length(pi0.grids))
        min.idx=which(starts.rslt==min(starts.rslt),TRUE)

    #cat("Done!",fill=TRUE)
        starts=c(pi0.grids[min.idx[1,2]], scale.grids[min.idx[1,1]])
    

    grad=function(parms){
        pi0=parms[1]; sc=parms[2]
        this.dt1=dt(tstat/sc,df)/sc
        ft=pi0*dt0+(1-pi0)*this.dt1
        d.dpi0=-2*sum((dt0-this.dt1)/ft)
        d.dsc= -2*sum((1-pi0)/sqrt(df)/beta(.5,df/2)*sc^(df-1)*(t2-sc*sc)*(t2/df+sc*sc)^(-df/2-1.5)/ft)
        c(d.dpi0,d.dsc)
    }


    #cat("Constrained Optimization...")        
    #    ui=matrix(c(1,0,
    #               -1,0,
    #                0,1),3,2,by=TRUE)
    #    ci=c(0,-1,1)
    #    (coptim.fit=constrOptim(starts,logLik,NULL,ui,ci))
    (optim.fit=optim(starts,logLik,grad,lower=c(1e-4, 1+1e-4),upper=c(1-1e-4, Inf),method='L-BFGS-B'))
    #cat("Done!",fill=TRUE)

    ans=optim.fit$par
    #    names(ans)=c('pi0','scale.fact')
    #    nparm=length(starts)
    #    attr(ans,'fit')=optim.fit
    ##    attr(ans,'dt1.delta')=dt1.delta
    #    attr(ans,'df')=df
    #    attr(ans,'equiv.sd.ncp')=sqrt(optim.fit$par[2]^2-1)
    #return(ans)


    lfdr.final=ans[1]*dt0
    dens=(lfdr.final+(1-ans[1])*dt(tstat/ans[2],df)/ans[2])
    lfdr.final=lfdr.final/dens
    final.scale=ans[2]
    log.dens=log(dens)
    crit=-log.dens+2/length(tstat)
    intercept=logit(1-ans[1])
    hess.nLogLik.pen=function(parms, return.dense=FALSE,...)  ## also depends on dt0, tstat, spar.Pen.mat, H
    {   r=parms[1]; scale.fact=1+exp(r);  # scale.fact=parms[1]; alternative parameterization to remove boundary
        betas=parms[-1]
        pi0s=drop(1/(1+exp(betas)))
        dt1=dt(tstat/scale.fact,df)/scale.fact
        f=drop(pi0s*dt0+(1-pi0s)*dt1)
        lfdrs=drop(pi0s*dt0/f)

        w.beta=drop((dt0-dt1)/f*(1-pi0s)*pi0s)
        W.beta.beta=w.beta*drop(pi0s*lfdrs-(1-pi0s)*(1-lfdrs))
        #        Wbb.H=W.beta.beta*H
        #        HWH=as(crossprod(H,Wbb.H),'symmetricMatrix')
        #        hess.bb=HWH+spar.Pen.mat
        hess.bb=sum(W.beta.beta)

        t2_s2=(tstat*tstat-scale.fact*scale.fact)
        t2vs2=(tstat*tstat/df+scale.fact*scale.fact)

        d.dr.i=drop(-t2_s2*exp(log(1-1/scale.fact)-log(t2vs2)+log(1-lfdrs)))
        w.br=drop(lfdrs*d.dr.i)
        d.dbr=sum(w.br)


        d2.dr2=sum((scale.fact-1)^2/scale.fact*(1-lfdrs)*t2_s2/t2vs2*(
                    -lfdrs/scale.fact*t2_s2/t2vs2+1/scale.fact/(1-scale.fact)
                    +2*scale.fact/t2_s2 +2*scale.fact/t2vs2
               ))
        
        
        #        ans=rBind(c(d2.dr2, drop(d.dbr)),cBind(drop(d.dbr),hess.bb))
        #        if(return.dense) as.matrix(ans) else as(ans, 'symmetricMatrix')
        matrix(c(d2.dr2, d.dbr, d.dbr, hess.bb), 2,2)
    }
    asym.vcov=solve(hess.nLogLik.pen(c(log(final.scale-1), logit(1-ans[1]))))

    ans=list(lfdr=lfdr.final,
             model=list(tstat=tstat, df=df, x=NULL), 
             scale.fact=list(scale.fact=final.scale, sd.ncp=sqrt(final.scale^2-1),  r=log(final.scale-1), 
                        t.cross=sqrt(df*(final.scale^(2/(df+1))-1)/(1-final.scale^(-2*df/(df+1))))),
             pi0=rep(ans[1], length(tstat)),
             tuning=list(mean=mean(crit),
                        var=var(log.dens)/length(tstat),  
                        grp=rep(1,length(tstat)),  method='NIC', 
                        final=mean(crit)),
             spar=list(all=Inf, final=Inf, final.idx=1),
             enp=list(raw=2, logistic=2, final=2, good.idx=TRUE),
             fit=list(intercept=intercept, 
                      covariate.idx=NULL,
                      f.covariate=NULL,
                      f=rep(intercept,length(tstat)),
                      beta=intercept,
                      H=matrix(1,length(tstat),1),
                      asym.vcov=asym.vcov),
             NPLL=list(NPLL=-sum(log.dens), logLik=log.dens, penalty=0, 
                       saturated.ll=#attr(logLik.saturated,'logLik'),
                                    ifelse(abs(tstat)<sqrt(df*(final.scale^(2/(df+1))-1)/(1-final.scale^(-2*df/(df+1)))), 
                                           dt(tstat,df,log=TRUE), dt(tstat/final.scale,df,log=TRUE)-log(final.scale))
                      )
            )

    class(ans)='hisemit'
    return(ans)


    if(FALSE){######## the following EM code is NOT corrected, but seems useless to use EM
        #        pi0.1=NA;parm.1=NA
        #        pi0.0=ans[1]
        #        parm.0=coptim.fit$par[[2]]
        #        logLik1=function(parms) sum((1-lfdr.0)*log(dt1.delta%*%
        #            (pncp(deltas[-1],parms)-pncp(deltas[-length(deltas)],parms))))
        #        K=proc.time()[3]
        #        repeat{
        ##            K=K+1
        #            f.0=pi0.0*dt0+(1-pi0.0)*dt1.delta%*%
        #                (pncp(deltas[-1],parm.0)-pncp(deltas[-length(deltas)],parm.0))
        #            lfdr.0=pi0.0*dt0/f.0
        #            parm.1=optimize(logLik1,2^c(-5,3),maximum=TRUE)$maximum
        #            if(FALSE){
        #                tmp=runif(50)
        #                tmpa=sapply(tmp,function(xx)-logLik(c(xx,parm.1)))
        #                tmp0=-logLik(c(pi0.0,parm.0))
        #                idx=tmpa>tmp0
        #                if(sum(idx>0)){cat("GEM",fill=TRUE);pi0.1=tmp[idx][which.max(abs(tmp[idx]-pi0.0))]}
        #                else pi0.1=mean(lfdr.0)
        #            }
        #            pi0.1=mean(lfdr.0)
        #            if(abs(parm.1-parm.0)<1e-5 && abs(pi0.1-pi0.0)<1e-5) break
        #            if(proc.time()[3]-K>1){
        #                cat('K=',proc.time()[3]-K,'\tpi0=',pi0.0,'\tgamma2=',parm.0,'\tlogLik=',-logLik(c(pi0.0,parm.0)),
        #                    '\tdelta pi0=',1e6*abs(pi0.0-pi0.1),'\tdelta gamma2=',1e6*abs(parm.0-parm.1),fill=TRUE)
        #                K=proc.time()[3]
        #            }                
        #            pi0.0=pi0.1
        #            parm.0=parm.1
        #        }
        #        ans.em=pi0.1
        ##        attr(ans.em,'gamma')=parm.1
        ##        attr(ans.em,'logLik')=-.5*logLik(c(pi0.1,parm.1))
        ##        return(list(ans,ans.em))
        #
        #        c(ans.em,parm.1)
    }
}



scaledTMix.psat=function(tstat,df, upper0=2)
{
    logLik=function(scale, takeSum=TRUE){
        n2ll=-2*
            pmax(dt(tstat,df,log=TRUE), dt(tstat/scale, df, log=TRUE)-log(scale))
      if(takeSum)sum(n2ll) else n2ll
    }
    if(upper0<=1)upper0=2
    upper=upper0
    lower=1
    repeat{
        optimize.rslt=optimize(logLik,c(lower,upper))
        if(optimize.rslt$minimum/upper<.99)break
        lower=1+(optimize.rslt$minimum-1)/upper0
        upper=upper*upper0
    }
    ans=optimize.rslt$minimum
    names(ans)='scale.fact'
    attr(ans,'equiv.sd.ncp')=sqrt(optimize.rslt$minimum^2-1)
    attr(ans,'df')=df
    attr(ans,'fit')=optimize.rslt
    attr(ans,'n2ll')=logLik(optimize.rslt$minimum,takeSum=FALSE)
    attr(ans,'pi0')=ifelse(dt(tstat,df)>dt(tstat/ans,df)/ans, 1, 0)
    ans
}

scaledTMix.sat=function(tstat,df)
{
    abs.t=abs(tstat)
    s=pmax(abs.t,1)
    pi0s=abs.t<=1
    logLik=ifelse(pi0s, dt(abs.t,df,log=TRUE), dt(1,df,log=TRUE)-log(abs.t))
    ans=s
    attr(ans,'pi0')=pi0s*1
    attr(ans,'logLik')=logLik

    ans
}


if(FALSE){### testing code



    G=25000
    sdncp=15.8
    n1=n2=10
    df=n1+n2-2
    x=1:G
    f=function(x)sin(x*pi/100)+1
    #Pi.i=rep(.7,G) 
    Pi.i=1/(1+exp(f(x)))
    mean(Pi.i)

    Z.i=rbinom(G,1,1-mean(Pi.i))    #rbinom(G,1,1-Pi.i)
    t0.i=rt(G,df)
    ncp.i=rnorm(G,0,sdncp)
    t1.i=rt(G,df,ncp.i)
    t.i=ifelse(Z.i==0,t0.i,t1.i)
    pvals=2*pt(abs(t.i),df,lower=FALSE)


    source("/Users/longor/MDS/scaledTMix.R")
    #debug(scaledTMix.null)
    (tmp.scl=scaledTMix.null(t.i,df))
    (tmp.sat=scaledTMix.sat(t.i,df))

    require(qvalue); qvalue(pvals)$pi0  
    require(ROCR); performance(prediction(abs(t.i),Z.i),'auc')@y.values[[1]]

    ######### below is not for this test
    #G.test=200
    #test.x=sort(runif(G.test,min(x),max(x)))
    #pi.test=1/(1+exp(f(test.x)))
    #Z.test=rbinom(G.test,1,1-pi.test)
    #t0.test=rt(G.test,df)
    #ncp.test=rnorm(G.test,0,sdncp)
    #t1.test=rt(G.test,df,ncp.test)
    #ti.test=ifelse(Z.test==0,t0.test,t1.test)
    #
    #
    #
    ##source("C:\\Users\\longor\\Desktop\\work\\MDS\\GOSim\\scaledTMix.R")
    ##source("Y:\\MDS\\scaledTMix.R")
    #source("/Users/longor/MDS/scaledTMix.R")
    #tmp=scaledTMix(t.i,df)
    ##tmp=penLik.aem(t.i,x,n1=n1, dt1.delta.all=dt1.delta)
    ##tmp=penLik.aem(t.i,x,n1=n1, dt1.delta.all=dt1.delta,cv.fold=10)
    ##tmp=penLik.aem(t.i,x,n1=n1,spar=1e3, dt1.delta.all=dt1.delta)
    #
    #
    ##
    ##summary(attr(tmp,'true.lfdr'))
    ##summary(attr(tmp,'lfdr'))
    ##cor(attr(tmp,'true.lfdr'),attr(tmp,'lfdr'))
    ##plot(attr(tmp,'true.lfdr'),attr(tmp,'lfdr'))
    #
    #t.pred=prediction(abs(t.i),Z.i)
    ##ideal.pred=prediction(-attr(tmp,'true.lfdr'),Z.i)
    #my.pred=prediction(-tmp[1:G],Z.i)
    #
    #performance(t.pred,'auc')@y.values[[1]]
    ##performance(ideal.pred,'auc')@y.values[[1]]
    #performance(my.pred,'auc')@y.values[[1]]
    #
    #plot(performance(t.pred,'tpr','fpr'),col=1)
    #plot(performance(ideal.pred,'tpr','fpr'),col=4,add=TRUE)
    #plot(performance(my.pred,'tpr','fpr'),col=2,add=TRUE)
    #
    #plot(attr(tmp,'true.lfdr'),attr(tmp,'lfdr')); abline(0,1,col=2,lwd=2)

}
