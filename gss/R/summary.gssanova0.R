## Summarize gssanova objects
summary.gssanova0 <- function(object,diagnostics=FALSE,...)
{
    y <- model.response(object$mf,"numeric")
    wt <- model.weights(object$mf)
    offset <- model.offset(object$mf)
    if ((object$family=="nbinomial")&(!is.null(object$nu))) y <- cbind(y,object$nu)
    dev.resid <- switch(object$family,
                        binomial=dev.resid.binomial(y,object$eta,wt),
                        nbinomial=dev.resid.nbinomial(y,object$eta,wt),
                        poisson=dev.resid.poisson(y,object$eta,wt),
                        inverse.gaussian=dev.resid.inverse.gaussian(y,object$eta,wt),
                        Gamma=dev.resid.Gamma(y,object$eta,wt),
                        weibull=dev.resid.weibull(y,object$eta,wt,object$nu),
                        lognorm=dev.resid.lognorm(y,object$eta,wt,object$nu),
                        loglogis=dev.resid.loglogis(y,object$eta,wt,object$nu))
    dev.null <- switch(object$family,
                       binomial=dev.null.binomial(y,wt,offset),
                       nbinomial=dev.null.nbinomial(y,wt,offset),
                       poisson=dev.null.poisson(y,wt,offset),
                       inverse.gaussian=dev.null.inverse.gaussian(y,wt,offset),
                       Gamma=dev.null.Gamma(y,wt,offset),
                       weibull=dev.null.weibull(y,wt,offset,object$nu),
                       lognorm=dev.null.lognorm(y,wt,offset,object$nu),
                       loglogis=dev.null.loglogis(y,wt,offset,object$nu))
    w <- object$w
    if (is.null(offset)) offset <- rep(0,length(object$eta))
    ## Residuals
    res <- 10^object$nlambda*object$c 
    ## Fitted values
    fitted <- object$eta
    fitted.off <- fitted-offset
    ## dispersion
    sigma2 <- object$varht
    ## RSS, deviance
    rss <- sum(res^2)
    dev <- sum(dev.resid)
    ## Penalty associated with the fit
    penalty <- sum(object$c*fitted.off*sqrt(w))
    penalty <- as.vector(10^object$nlambda*penalty)
    ## Calculate the diagnostics
    if (diagnostics) {
        ## Obtain retrospective linear model
        comp <- NULL
        for (label in object$terms$labels) {
            if (label=="1") next
            if (label=="offset") next
            comp <- cbind(comp,predict(object,object$mf,inc=label))
        }
        comp <- cbind(comp,yhat=fitted.off,y=fitted.off+res/sqrt(w),e=res/sqrt(w))
        term.label <- object$terms$labels[object$terms$labels!="1"]
        term.label <- term.label[term.label!="offset"]
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
        if (qr(corr)$rank<dim(corr)[2]) kappa <- rep(Inf,len=dim(corr)[2])
        else kappa <- as.numeric(sqrt(diag(solve(corr))))
        ## Obtain decomposition of penalty
        rough <- as.vector(10^object$nlambda*t(comp[,term.label])%*%object$c/penalty)
        names(kappa) <- names(rough) <- term.label
    }
    else decom <- kappa <- cosines <- rough <- NULL
    ## Return the summaries
    z <- list(call=object$call,family=object$family,method=object$method,iter=object$iter,
              fitted=fitted,dispersion=sigma2,residuals=res/sqrt(w),rss=rss,
              deviance=dev,dev.resid=sqrt(dev.resid)*sign(res),nu=object$nu,
              dev.null=dev.null,penalty=penalty,
              pi=decom,kappa=kappa,cosines=cosines,roughness=rough)
    class(z) <- "summary.gssanova0"
    z
}
