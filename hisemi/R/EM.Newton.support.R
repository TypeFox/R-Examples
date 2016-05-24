plot.hisemit=function(x, type=c('tuning', 'residual'), ...)
{
    type=match.arg(type)
    if(type=='tuning') {
        plotHisemitTuning(x,...)
    }else if(type=='residual'){
        plotHisemitResid(x,...)
    }else {
        warning("plot type not implemented")
    }
}
print.hisemit=function(x,...)
{
#    print(unclass(x),...)
    cat('sd.ncp =', x$scale.fact$sd.ncp, fill=TRUE)
    cat("enp =", x$enp$final, fill=TRUE)
    cat("smoothing parameter =", x$spar$final, fill=TRUE)
    cat("logLik =", sum(x$NPLL$logLik), ',\tpenalty =',x$NPLL$penalty, fill=TRUE)
    cat('lfdr summary: ',fill=TRUE)
    print(summary(x$lfdr))
    cat('pi0 summary: ', fill=TRUE)
    print(summary(x$pi0))

    invisible(x)

}

summary.hisemit=function(object,...)
{
    print(object,...)
}

plotHisemitResid=function(obj, y.type=c('hist','scatter'), x.type=c('lfdr','pi0','f'),...)
{
    y.type=match.arg(y.type)
    if(y.type=='hist')        return(hist(resid(obj),xlab='residuals',main='Histogram of Residuals',...))

    if(y.type=='scatter'){
        x.type=match.arg(x.type)
        x.fit=fitted(obj, fitted.type=x.type)
        return ( plot(x.fit, resid(obj,...), xlab=x.type, ylab='residual', pch='.') )
    }
    {
        stop("unimplenmented residual plot") # shouldn't reach this line
    }
}

plotHisemitTuning=function(obj, SE=FALSE, add=FALSE,  ...)
{   if(length(setdiff(obj$spar$all,Inf))<1) stop('need finite smoothing parameters')
    tuning.method=obj$tuning$method
    spar.exp=log10(obj$spar$all)
    final.spar.exp=log10(obj$spar$final)
    criterion.mean=obj$tuning$mean
    enps=colMeans(obj$enp$logistic)
    enps0=colMeans(obj$enp$raw)

    criterion.var=obj$tuning$var
    wt.cv=table(obj$tuning$grp)
    G=length(obj$lfdr)
    criterion.mean.mean=crossprod(wt.cv, criterion.mean)/sum(wt.cv)
    criterion.var.se=sqrt(crossprod(wt.cv*wt.cv, criterion.var)/G/G)


    goodenp.idx=obj$enp$good.idx
    s.mode.i=which(goodenp.idx)[1]
    imin.cv=s.mode.i-1+which.min(criterion.mean.mean[goodenp.idx])


    myplot=if(add)lines else plot

    crit.range=range(criterion.mean.mean, na.rm=TRUE)*G
    enp.range=c(2+ncol(obj$model$x)*(obj$model$pen.order-1), max(obj$enp$logistic))
            
    if(!add){    par(mar=c(5,4,4,4)+.1)}

    if(SE){
        myplot(spar.exp[goodenp.idx], G*drop(criterion.mean.mean)[goodenp.idx], ylab=if(add)''else tuning.method, 
            xlab='log10(smoothing par.)',
            ylim=range(G*c(criterion.mean.mean+criterion.var.se,criterion.mean.mean-criterion.var.se)),
                type='o',...)
        for(i in 1:length(spar.exp))
            lines(c(spar.exp[i],spar.exp[i]),G*(criterion.mean.mean[i]+c(1,-1)*criterion.var.se[i]),col=2,lwd=3)
        abline(v=spar.exp[imin.cv])
    #        spar.exp[criterion.mean.mean<(criterion.mean.mean[imin.cv]+criterion.var.se[imin.cv])]=Inf
    #        cv.1se=tail(which.max(spar.exp),1)
    #        abline(v=spar.exp[cv.1se],lty=2)
    #        axis(1,spar.exp[cv.1se])
    #    abline(h=criterion.mean.mean[imin.cv]+criterion.var.se[imin.cv])
    }else{

            myplot(spar.exp,G*drop(criterion.mean.mean), ylab=if(add)''else tuning.method, xlab='log10(smoothing par.)',
                 type='o', ...)
            abline(v=final.spar.exp)

            lines(spar.exp, (enps-mean(enp.range))/diff(enp.range)*diff(crit.range)+mean(crit.range), col=4, lwd=2)
            points(spar.exp, (enps0-mean(enp.range))/diff(enp.range)*diff(crit.range)+mean(crit.range), pch=3, 
                    col=ifelse(goodenp.idx, 'blue', 'gray'))
            right.tick=round(c(enps[which(spar.exp==final.spar.exp)], max(enps,na.rm=TRUE),
                            if(!is.na(enps[1])) enps[1] else enps[s.mode.i-1], enp.range))
            axis(4, (right.tick-mean(enp.range))/diff(enp.range)*diff(crit.range)+mean(crit.range), right.tick, 
                    col=4, col.ticks=4, lwd=2, col.axis=4)
            mtext('eff. #parms', 4, col=4,line=2)
            abline(v=spar.exp[s.mode.i],lty=2)

    }
    axis(1,spar.exp[imin.cv])
    axis(1,range(spar.exp))
    
    invisible(crit.range)
}


coef.hisemit=function(object, scale.parameterization=c('r','scale.factor','sd.ncp'), ...)
{
    scale.parameterization=match.arg(scale.parameterization)
    ans=c(if(scale.parameterization=='r') {
                object$scale.fact$r
          }else if (scale.parameterization=='scale.factor') {
                object$scale.fact$scale.fact
          }else if (scale.parameterization=='sd.ncp') {
                object$scale.fact$sd.ncp
          }else {NA}, 
          object$fit$beta)
    names(ans)=c(scale.parameterization, paste("beta", 0:(length(ans)-2),sep='.'))
    ans
}

fitted.hisemit=function(object, fitted.type=c('lfdr','fpp','pi0','f'), gene.list, component, ...)
{
    fitted.type=match.arg(fitted.type)
    if(fitted.type=='lfdr')        return(object$lfdr)
    if(fitted.type=='fpp') {
        l=object$lfdr
        ord.l=order(l)
        r=1-logit.inv(cumsum(logit(1-l[ord.l]))/(seq(along=l)))
        ans=l
        ans[ord.l]=r
        return (ans)
    }
    if (fitted.type=='pi0')         return(object$pi0)
    if (fitted.type=='f') {
        if(missing(component))  return(object$fit$f)
        if(!is.na(pmatch(component,'intercept')) || (length(component)==1 && component==0)) { 
            return(rep(object$fit$intercept, length(object$lfdr)))
        }
        if (is.numeric(component)) {
            stopifnot(all(round(component)<=ncol(object$fit$f.covariate)))
            component=sort(unique(  pmax(pmin(round(component), ncol(object$fit$f.covariate)),0) ))
            if(length(component)==0) stop("invalid component")
            if(component[1]==0){ 
                ans=matrix(object$fit$intercept, length(object$lfdr), ncol=1, dimnames=list(NULL, 'intercept'))
                component=component[-1]
                if(length(component)==0) return(drop(ans))
            }else ans=matrix(NA_real_, nrow(object$fit$f.covariate),0)
            ans.rest=object$fit$f.covariate[,component,drop=FALSE]
            colnames(ans.rest)=paste("covariate",component,sep='.')
            return(drop(cbind(ans, ans.rest)))
        }
        {
            warning("component argument is unknown. overall fit is returned")
            return(object$fit$f)
        }
    }
    {
        stop("the program shouldn't reach this line")
    }
}

residuals.hisemit=function(object, residual.type='deviance', ...)
{
    residual.type=match.arg(residual.type)
    if(residual.type=='deviance'){
        return (sign(object$model$tstat)*sqrt(2*(object$NPLL$saturated.ll-object$NPLL$logLik))  )
    }
    {
        stop("residuals other than deviance type are not implemented")
    }

}

logLik.hisemit=function(object, take.sum=TRUE, ...)
{
    ans=object$NPLL$logLik
    enp=object$enp$final

    if(take.sum) ans=sum(ans)
    attr(ans, 'df')=enp
    class(ans)='logLik'
    ans
}

logit=make.link("logit")$linkfun
logit.inv=make.link("logit")$linkinv

vcov.hisemit=function(object,...)
{
    object$fit$asym.vcov
}

#NIC=function(obj,...)
#{
#    AIC(obj,...)
#}

confint.hisemit=function(object, parm=c('lfdr', 'fpp', 'beta', 'scale.fact','sd.ncp','r','coef','pi0','f'),level=.95,component,... )
{
    parm=match.arg(parm)
    tstat=object$model$tstat
    df.resid=length(tstat)-object$enp$final
    t.1malpha=qt(0.5+level/2, df.resid)
    V=vcov(object)
    s=object$scale.fact$scale.fact
    df=object$model$df
    if(parm=='lfdr')    {
        l=fitted(object,'lfdr')
        d.logit.1ml.dr=(tstat*tstat-s*s)*df*(s-1)/(tstat*tstat+df*s*s)/s
        grad=cBind(d.logit.1ml.dr,object$fit$H)

        ########### this is toooooooooo slow
        #        var.logit.1ml=sapply(1:length(tstat), function(i)drop(grad[i,,drop=FALSE]%*%V%*%grad[i,]))
        ########### this is much faster
        chol.V=chol(nearPD(V)$mat)  ## upper triagular
        gradL=tcrossprod(grad,chol.V)
        var.logit.1ml=rowSums(gradL*gradL)
        ########### END: this is much faster

        CL.logit.1ml=logit(1-l)+outer( t.1malpha*sqrt(var.logit.1ml), c(1,-1) )
        CL.l=1-logit.inv(CL.logit.1ml)
        colnames(CL.l)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        rownames(CL.l)=names(tstat)
        return(CL.l)
    }
    if(parm=='beta')    {
        var.beta=diag(V)[-1]
        cl.beta=coef(object)[-1]+outer(t.1malpha*sqrt(var.beta), c(-1,1))
        colnames(cl.beta)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        rownames(cl.beta)=names(coef(object))[-1]
        return(cl.beta)
    }
    if(parm=='scale.fact'){
        var.r=drop(V[1,1])
        cl.r=object$scale.fact$r+c(-1,1)*t.1malpha*sqrt(var.r)
        cl.s=1+exp(cl.r)
        names(cl.s)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        return(cl.s)
    }
    if(parm=='sd.ncp'){
        var.r=drop(V[1,1])
        cl.r=object$scale.fact$r+c(-1,1)*t.1malpha*sqrt(var.r)
        cl.sdncp=sqrt((1+exp(cl.r))^2-1)
        names(cl.sdncp)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        return(cl.sdncp)
    }
    if(parm=='r'){
        var.r=drop(V[1,1])
        cl.r=object$scale.fact$r+c(-1,1)*t.1malpha*sqrt(var.r)
        names(cl.r)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        return(cl.r)
    }
    if(parm=='coef'){
        var.coef=diag(V)
        cl.coef=coef(object)+outer(t.1malpha*sqrt(var.coef), c(-1,1))
        colnames(cl.coef)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        rownames(cl.coef)=names(coef(object))
        return(cl.coef)
    }
    if(parm=='f'){
        if(missing(component)){
            cov.idx=seq_len(ncol(object$fit$H))
            fitted.f=fitted(object, 'f')
        }else {
            if(!is.na(pmatch(component,'intercept')) || (length(component)==1 && component==0)) { 
                return(rep(1,length(object$model$tstat))%o%confint(object, 'beta')[1,])
            }
            if (is.numeric(component)) {
                component=sort(unique(  pmax(pmin(round(component), ncol(object$fit$f.covariate)),0) ))
                if(length(component)!=1) stop("invalid component")
                if(component[1]==0){ return(rep(1,length(object$model$tstat))%o%confint(object, 'beta')[1,]) }
                cov.idx=which(object$fit$covariate.idx==component)# this does not count the first scale parameter
                fitted.f=fitted(object, 'f', component=component)
            }
        }
        chol.V=chol(nearPD(V)$mat)
        grad=cBind(0, object$fit$H)
        grad[,-cov.idx-1L]=0
        gradL=tcrossprod(grad,chol.V)
        var.f=rowSums(gradL*gradL)
        CL.f=fitted.f+outer( t.1malpha*sqrt(var.f), c(-1,1) )
        colnames(CL.f)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        rownames(CL.f)=names(tstat)
        return(CL.f)
    }
    if(parm=='pi0'){
        chol.V=chol(nearPD(V)$mat)
        grad=cBind(0, object$fit$H)
        gradL=tcrossprod(grad,chol.V)
        var.f=rowSums(gradL*gradL)

        CL.f=fitted(object,'f')+outer( t.1malpha*sqrt(var.f), c(1,-1) )
        CL.pi0=1-logit.inv(CL.f)
        colnames(CL.pi0)=paste(round(c(1-0.5-level/2, 0.5+level/2)*100, 3), '%',sep='')
        rownames(CL.pi0)=names(tstat)
        return(CL.pi0)
    }
    if(parm=='fpp') {
        stop('not yet implemented')
    }

    {
        stop("shouldn't reach this line")
    }

}

