EMupdate=function(starts, 
    nLogLik.pen,
    optim.method,
    H,
    tstat,
    df,
    dt0,
    spar.Pen.mat,
    em.iter.max=10,
    em.beta.iter.max=1,
    scale.conv=1e-3,
    lfdr.conv=1e-3,
    NPLL.conv=1e-3,
    debugging=FALSE
){
    if(em.iter.max<1) return(starts)
    
    sf.start=1+exp(starts[1])
    beta.start=starts[-1]


    logLik.scale=function(scale.fact) sum((1-lfdr.old)*(dt(tstat/scale.fact,df,log=TRUE)-log(scale.fact)))
    max.tstat=max(abs(tstat))
    quadsol.logLik.scale=function(mle=FALSE,debugging=FALSE,...){    ### optimization function for scale.factor
        if(!mle){ 
            ds=function(s) sum( (1-lfdr.old)*(tstat^2-s^2)*df/((tstat^2+ df* s^2 )*s))
            optim.rslt=optim(scale.old, logLik.scale, ds, method='L-BFGS-B', lower=1, upper=max.tstat, 
                control=list(fnscale=-1,maxit=em.beta.iter.max))
            return(optim.rslt$par)
        }
        if(debugging) cat("optimizing wrt scale factor...", fill=TRUE)
        lower=1
        upper0=2
        upper=upper0
        repeat{
            optimize.rslt=optimize(logLik.scale,c(lower,upper),maximum=TRUE,tol=scale.conv)
            if(optimize.rslt$maximum/upper<.99)break
            lower=1+(optimize.rslt$maximum-1)/upper0
            upper=upper*upper0
        }
        optimize.rslt$maximum
    }

    em.nr.obj=function(betas, with.grad.hess=FALSE){
        fx=drop(H%*%betas)
        pi0s=1/(1+exp(fx))
        ans=sum(-(1-lfdr.old)*fx-log(pi0s))+.5*drop(betas%*%spar.Pen.mat%*%betas)
        if(with.grad.hess){
            attr(ans,'gradient')=drop((-pi0s+lfdr.old)%*%H + betas%*%spar.Pen.mat)
            w.h.tilde=pi0s*(1-pi0s)
            WH=H*sqrt(w.h.tilde)                              ## this is equivalent to diag(wt)%*%H but faster
            HWH=crossprod(WH)                    ### this crossprod is the slowest part
            attr(ans,'hessian')=as.matrix(HWH + spar.Pen.mat)
        }
        ans
    }
    em.nr.grad=function(betas){
        pi0s=drop(1/(1+exp(H%*%betas)))
        drop((-pi0s+lfdr.old)%*%H + betas%*%spar.Pen.mat)
    }
    em.nr.hess=function(betas){
        pi0s=drop(1/(1+exp(H%*%betas)))
        w.h.tilde=pi0s*(1-pi0s)
        WH=H*sqrt(w.h.tilde)                              ## this is equivalent to diag(wt)%*%H but faster
        HWH=crossprod(WH)                    ### this crossprod is the slowest part
        as.matrix(HWH + spar.Pen.mat)
    }


################ EM algorithm loop
    ###### starting values for EM
    beta.old=beta.start
    scale.old=sf.start
    fx.old=drop(H%*%beta.old)  ## vector
    pi0.old=1/(1+exp(fx.old))
        dt1=dt(tstat/scale.old,df)/scale.old
        lfdr.old <- pi0.old*dt0
    f.old=lfdr.old+(1-pi0.old)*dt1
    lfdr.old=lfdr.old/f.old
    npll.old=nLogLik.pen(c(log(scale.old-1),beta.old))
if(debugging) cat('\n',0,'npll=',npll.old,'scale.old=',scale.old, fill=TRUE)
if(debugging) {scale.time=0; newton.time=0}
    em.iter=1
    while(em.iter<=em.iter.max){
################# this is approximately optimizing scale.fact
        if(debugging) t0=proc.time()[3]
        scale.new=quadsol.logLik.scale(FALSE,debugging)
        if(debugging) scale.time=scale.time+proc.time()[3]-t0

################# this is optimizing  beta using Newton-type alagorithms
        if(debugging) t0=proc.time()[3]
        if(optim.method%in%c("BFGS","CG","L-BFGS-B","Nelder-Mead", "SANN")){ ## call optim
            optim.fit=optim(beta.old, em.nr.obj, em.nr.grad, method=optim.method, 
                        control=list(maxit=em.beta.iter.max, trace=if(debugging)3 else 0, REPORT=30))
            beta.new=optim.fit$par
        }else if (optim.method=='nlminb') {
            nlminb.fit=nlminb(beta.old, em.nr.obj, em.nr.grad, em.nr.hess, #return.dense=TRUE, 
                        control=list(eval.max=2000,iter.max=em.beta.iter.max))
            beta.new=nlminb.fit$par
        }else if (optim.method=='NR'){
            NR.fit=NRupdate(em.nr.obj, beta.old, iter.max=em.beta.iter.max, with.grad.hess=TRUE)
            beta.new=NR.fit
        }
        if(debugging) newton.time=newton.time+proc.time()[3]-t0


        fx.new=drop(H%*%beta.new)  ## vector
        pi0.new=1/(1+exp(fx.new))
        dt1=dt(tstat/scale.new, df)/scale.new
        lfdr.new <- pi0.new*dt0
        f.new=lfdr.new+(1-pi0.new)*dt1
        lfdr.new=lfdr.new/f.new
        npll.new=nLogLik.pen(c(log(scale.new-1),beta.new))


if(debugging) cat('\n',em.iter,'npll=',npll.new,'\tdiff.npll',npll.old-npll.new, '\tscale.new=',scale.new, 
    '\tdiff.scale=',(scale.new-scale.old),'\tdiff.li=',max(abs(lfdr.new-lfdr.old)), fill=TRUE )
################ convergence checking
        if( abs(scale.new-scale.old)<scale.conv && 
            max(abs(lfdr.old-lfdr.new))<lfdr.conv &&
            (npll.old-npll.new)<NPLL.conv) break 

###############  setting conditions for next iteration
        pi0.old<-pi0.new    ## not necessary
        scale.old=scale.new
        beta.old=beta.new
        lfdr.old=lfdr.new
        f.old=f.new         ## not necessary
        npll.old=npll.new
        fx.old=fx.new       ## not necessary

   em.iter=em.iter+1
   } ## of EM loop
   if(debugging)cat("\tscaling.time=", scale.time, "\tnewton.EM.time=", newton.time, fill=TRUE)

   ans=c(log(scale.new-1), beta.new)
   ans
}
