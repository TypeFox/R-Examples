# fit.4pl=FALSE; var.model="power"; robust="mean"; method="drm"; max.iter=20; reltol=1e-3; verbose=TRUE; log.both.sides=F
glspl=function (formula, data, max.iter=50, reltol=1e-3, fit.4pl=FALSE, gof.threshold=0.2, verbose=FALSE) {

    # reorder, b/c drm 5pl is sensitive to the ordering of rows
    predictor.coln=all.vars(formula)[2]    
    data=data[order(data[[predictor.coln]]),] 
    
    outcome.var=model.frame(formula, data)[[1]]
    data[["outcome.var"]]=outcome.var # this is necessary, b/c drc uses model.weights to extract weight from a model.frame

    # set starting values
    tmpfit=drm.fit(formula, data=data, fit.4pl=fit.4pl, w=outcome.var^(-1), verbose=verbose>=2, gof.threshold=gof.threshold) # use model.frame which can take care of log transformation of outcome
    startVal = coef(tmpfit)
    startVal=cla2gh(startVal)
    
    fitte=fitted(tmpfit)
    resi = outcome.var-fitte 
    sigmasq.hat= getVarComponent.drc(tmpfit)
    g.hat=-1.2 * -2
    dev.0=sum( g.hat*log(fitte) + log(sigmasq.hat) + resi^2/(fitte^g.hat*sigmasq.hat) )
    old.fit=tmpfit
    if (verbose) cat("Iter 0. theta:", startVal, "g.hat:", g.hat, " deviance:", dev.0, "\n")
    
    # iterate between drm fit given a fixed weights and estimating the variance component parameters
    iterations=0
    theta=startVal        
    devs=c()
    while (TRUE) {
        iterations=iterations+1
        
        # estimate variance component, g.hat and sigmasq.hat, via uniroot
        f1=function(x,c,d,g,h,f) (c + (d - c)/(1 + exp(-log(f) - h/(d - c) * (1 + 1/f)^(1 + f) * (log(x) - g)))^f)
        fitte=f1(data[[predictor.coln]], theta["c"],theta["d"],theta["g"],theta["h"],ifelse(fit.4pl,1,theta["f"]))
        resi = outcome.var-fitte
        
        g.hat = try(uniroot(f=function(g,...) sum( 
            log(fitte) * ( 1 - resi^2/{fitte^g * mean(resi^2/fitte^g)} )
        ), interval=c(-3,3), fitte)$root, silent=TRUE) 
        # try one more time if an error happens, this time widen the interval
        if (inherits(g.hat, "try-error")) {
            g.hat = try(uniroot(f=function(g,...) sum( 
                log(fitte) * ( 1 - resi^2/{fitte^g * mean(resi^2/fitte^g)} ) # mean(resi^2/fitte^g) is sigmasq.hat
            ), interval=c(-10,10), fitte)$root) 
            if (inherits(g.hat, "try-error")) {
                cat("fail to find gamma, return null\n")
                return (NULL)
            }
        }

        w=fitte^-(g.hat/2) # /2 is needed b/c weight will be squared
        sigmasq.hat = mean(resi^2/fitte^g.hat)
        
        # estimate curve fit via drm
        fit= drm.fit(formula=formula, data=data, fit.4pl=fit.4pl, w=w, verbose=verbose>=2, gof.threshold=gof.threshold) 
        new.theta=cla2gh(coef(fit))
        # associate g.hat and sigmasq.hat with fit
        fit$var.power = g.hat
        fit$sigmasq.hat = sigmasq.hat
        
        # compute deviance
        fitte=fitted(fit)
        resi = outcome.var-fitte 
        devs = c(devs, sum( g.hat*log(fitte) + log(sigmasq.hat) + resi^2/(fitte^g.hat*sigmasq.hat) )  )
        
        if (verbose) cat("Iter "%+%iterations%+%".", "theta:", new.theta, "g.hat:", g.hat, " deviance:", devs[iterations], "\n")
        
        # stopping rules
        if (max(abs(1 - new.theta/theta)) < reltol) {
            if (verbose) cat("converged\n")
            fit$converged=TRUE
            break;
        } else if (iterations>=max.iter) {
            if (verbose) cat("max iter reached\n")
            fit$converged=FALSE
            break;
        } else if (devs[iterations]>ifelse(iterations>1,devs[iterations-1],dev.0)) {
            if (verbose) cat("dev increasing, stop and revert to previous fit\n")
            fit$converged=FALSE
            fit=old.fit 
            break;
        } else {
            # update theta, save fit to old.fit, and continue
            theta = new.theta 
            old.fit=fit
        }        
        
    } # end loop
    
    fit$iter=iterations
    
#        # try to see if can converge by gnls. it seems to lead to increased deviance
#        if (!fit$converged) {
#            fit.1 = gnls.fit (formula, data, fit.4pl=fit.4pl, verbose=verbose, startVal=coef(fit), varFun=varPower(g.hat/2))
#            if (!is.null(fit.1)) {
#                theta=coef(fit.1)
#                fitte.1=f1(data[[outcome.coln]], data[[predictor.coln]], theta["c"],theta["d"],theta["g"],theta["h"],theta["f"])
#                g.1=as.numeric(fit.1$modelStruct)*2
#                dev.1 = sum( g.1*log(fitte.1) + (data[[outcome.coln]]-fitte.1)^2/fitte.1^g.1 )
#
#            }
#        }

    fit
    
}
