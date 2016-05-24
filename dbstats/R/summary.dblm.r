

######################
#### summary.dblm ####
######################

 ## Description:
 ##     Summary of an object of class dblm.
 ##         - calculate the R-squared
 ##         - the adjusted R-squared
 ##         - the weighted residuals
 ##        


summary.dblm <-function(object,...)
{
    # recover attributes rdf, weights, residuals of dblm
    z <- object
    p <- z$eff.rank
    rdf <- z$df.residuals
    w <- z$weights/sum(z$weights)
    r <- z$residuals
    y0 <- z$y - sum(w*z$y)

    # weighted residuals
    r <- sqrt(z$weights) * r
    
    # residual sigma
    sigma<-attr(z,"sigma")
    
    # gvect
    gvec<-diag(attr(object,"G"))
    
    # Gemoetric variability
    gvar=t(w)%*%as.matrix(gvec)
    
    # Rsquared and Adj.RSquared
    R2 <- 1 - sum(w*(z$fitted.values-z$y)^2)/(sum(w*(y0^2)))
    adjR2 <- 1-  sum(w*(z$fitted.values-z$y)^2)/(sum(w*(y0^2)))*(length(z$y)-1)/(length(z$y)-1-p)
    rdf<-z$rdf
    
    # matched call
    call <-object$call
   
    crit.value<-NULL    
    if (attr(object,"method")=="OCV")
       crit.value<-object$ocv
    if (attr(object,"method")=="GCV")
       crit.value<-object$gcv
    if (attr(object,"method")=="AIC")
       crit.value<-object$aic
    if (attr(object,"method")=="BIC")
       crit.value<-object$bic

   
    
    # list to be returned
    ans <- list(residuals=r,sigma=sigma,r.squared=R2,adj.r.squared=adjR2,
            df.residuals=rdf,call=z$call,gvar=gvar,gvec=gvec,
            method=attr(object,"method"),eff.rank=object$eff.rank,
            rel.gvar=object$rel.gvar,crit.value=crit.value)
    class(ans)<-"summary.dblm"
    
    return(ans)
}
    
    
    
    