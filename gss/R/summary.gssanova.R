## Summarize gssanova objects
summary.gssanova <- function(object,diagnostics=FALSE,...)
{
    y <- model.response(object$mf,"numeric")
    wt <- model.weights(object$mf)
    offset <- model.offset(object$mf)
    if ((object$family=="nbinomial")&(!is.null(object$nu))) y <- cbind(y,object$nu)
    dev.null <- switch(object$family,
                       binomial=dev.null.binomial(y,wt,offset),
                       nbinomial=dev.null.nbinomial(y,wt,offset),
                       poisson=dev.null.poisson(y,wt,offset),
                       Gamma=dev.null.Gamma(y,wt,offset),
                       weibull=dev.null.weibull(y,wt,offset,object$nu),
                       lognorm=dev.null.lognorm(y,wt,offset,object$nu),
                       loglogis=dev.null.loglogis(y,wt,offset,object$nu))
    w <- object$w
    if (is.null(offset)) offset <- rep(0,length(object$eta))
    ## Residuals
    res <- residuals(object)*sqrt(w)
    dev.resid <- residuals(object,"deviance")
    ## Fitted values
    fitted <- fitted(object)
    ## dispersion
    sigma2 <- object$varht
    ## RSS, deviance
    rss <- sum(res^2)
    dev <- sum(dev.resid^2)
    ## Penalty associated with the fit
    obj.wk <- object
    obj.wk$d[] <- 0
    if (!is.null(model.offset(obj.wk$mf))) obj.wk$mf[,"(offset)"] <- 0
    penalty <- sum(obj.wk$c*predict(obj.wk,obj.wk$mf[object$id.basis,]))
    penalty <- as.vector(10^object$nlambda*penalty)
    if (!is.null(object$random)) {
        p.ran <- t(object$b)%*%object$random$sigma$fun(object$zeta,object$random$sigma$env)%*%object$b
        penalty <- penalty + p.ran
    }
    ## Calculate the diagnostics
    if (diagnostics) {
        ## Obtain retrospective linear model
        comp <- NULL
        p.dec <- NULL
        for (label in object$terms$labels) {
            if (label=="1") next
            if (label=="offset") next
            comp <- cbind(comp,predict(object,object$mf,inc=label))
            jk <- sum(obj.wk$c*predict(obj.wk,obj.wk$mf[object$id.basis,],inc=label))
            p.dec <- c(p.dec,10^object$nlambda*jk)
        }
        term.label <- object$terms$labels[object$terms$labels!="1"]
        term.label <- term.label[term.label!="offset"]
        if (!is.null(object$random)) {
            mf <- object$mf
            mf$random <- object$random$z
            comp <- cbind(comp,predict(object,mf,inc=NULL))
            p.dec <- c(p.dec,p.ran)
            term.label <- c(term.label,"random")
        }
        fitted.off <- fitted-offset
        comp <- cbind(comp,yhat=fitted.off,y=fitted.off+res/sqrt(w),e=res/sqrt(w))
        if (any(outer(term.label,c("yhat","y","e"),"==")))
            warning("gss warning in summary.gssanova: avoid using yhat, y, or e as variable names")
        colnames(comp) <- c(term.label,"yhat","y","e")
        ## Sweep out constant
        comp <- sqrt(w)*comp - outer(sqrt(w),apply(w*comp,2,sum))/sum(w)
        ## Obtain pi
        comp1 <- comp[,c(term.label,"yhat")]
        decom <- t(comp1) %*% comp1[,"yhat"]
        names(decom) <- c(term.label,"yhat")
        decom <- decom[term.label]/decom["yhat"]
        ## Obtain kappa, norm, and cosines        
        corr <- t(comp)%*%comp
        corr <- t(corr/sqrt(diag(corr)))/sqrt(diag(corr))
        norm <- apply(comp,2,function(x){sqrt(sum(x^2))})
        cosines <- rbind(corr[c("y","e"),],norm)
        rownames(cosines) <- c("cos.y","cos.e","norm")
        corr <- corr[term.label,term.label,drop=FALSE]
        if (qr(corr)$rank<dim(corr)[2])
            kappa <- rep(Inf,len=dim(corr)[2])
        else kappa <- as.numeric(sqrt(diag(solve(corr))))
        ## Obtain decomposition of penalty
        rough <- p.dec / penalty
        names(kappa) <- names(rough) <- term.label
    }
    else decom <- kappa <- cosines <- rough <- NULL
    ## Return the summaries
    z <- list(call=object$call,family=object$family,alpha=object$alpha,
              fitted=fitted,dispersion=sigma2,residuals=res/sqrt(w),rss=rss,
              deviance=dev,dev.resid=dev.resid,nu=object$nu,
              dev.null=dev.null,penalty=penalty,
              pi=decom,kappa=kappa,cosines=cosines,roughness=rough)
    class(z) <- "summary.gssanova"
    z
}
