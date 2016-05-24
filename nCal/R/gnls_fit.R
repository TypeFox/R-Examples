if(getRversion() >= "2.15.1")  utils::globalVariables(c("form.gnls"))

gnls.fit=function (formula, data, fit.4pl=FALSE, startVal=NULL, varFun=nlme::varPower(), verbose=FALSE) {
    
    outcome.coln=all.vars(formula)[1]
    predictor.coln=all.vars(formula)[2]    
    data=data[order(data[[predictor.coln]]),] # order the rows
    
    if (is.null(startVal)) {
        tmpfit=drm.fit(formula, data=data, fit.4pl=fit.4pl, w=(data[[outcome.coln]])^(-1), verbose=verbose)
        startVal = coef(tmpfit)
        startVal=cla2gh(startVal)
        if(verbose) print(startVal)        
    }
    if (!"h" %in% names(startVal)) startVal=cla2gh(startVal)
    
    control=nlme::gnlsControl(
        nlsTol=1e-1,  # nlsTol seems to be important, if set to 0.01, then often does not converge
        tolerance=1e-4, 
#        msTol=1e-1, 
#        minScale=1e-1, 
#        .relStep=1e-7,
        returnObject=TRUE, # allow the return of the fit when the max iter is reached
        maxIter=500, nlsMaxIter=50, opt="nlminb", optimMethod="BFGS", msVerbose=verbose>=2)

    # g-h parameterization works better than classical in convergence for some datasets, but still has difficulty with a lot of datasets
    if (!fit.4pl){   
        #eval(eval(substitute(expression( form.gnls <<- as.formula(outcome.coln%+%" ~ c + (d - c)/(1 + exp(b * (log("%+%predictor.coln%+%") - loge)))^f") ))))              
        #eval(eval(substitute(expression( form.gnls <<- as.formula(outcome.coln%+%" ~ exp(logc) + (exp(logd) - exp(logc))/(1 + exp(-log(f)-h/(exp(logd)-exp(logc))*(1+1/f)^(1+f)*(log("%+%predictor.coln%+%")-g)))^f") ))))              
        eval(eval(substitute(expression( form.gnls <<- as.formula(outcome.coln%+%" ~ c + (d - c)/(1 + exp(-log(f)-h/(d-c)*(1+1/f)^(1+f)*(log("%+%predictor.coln%+%")-g)))^f") ))))              
    } else {
        #eval(eval(substitute(expression( form.gnls <<- as.formula(outcome.coln%+%" ~ c + (d - c)/(1 + exp(b * (log("%+%predictor.coln%+%") - loge)))") ))))     
        eval(eval(substitute(expression( form.gnls <<- as.formula(outcome.coln%+%" ~ c + (d - c)/(1 + exp(-h/(d-c)*4*(log("%+%predictor.coln%+%")-g)))") ))))              
    }        
    
    fit= try(nlme::gnls(form.gnls, data=data, start=startVal, control=control, weights=varFun, na.action=na.omit), silent=FALSE)       
    if (class(fit)[1]=="try-error") {
        if (verbose) cat("return null fit\n")
        fit = NULL
    } else {
        vcov.=vcov(fit)
        if(is.null(vcov.)) {
            bad.se = T 
        } else if(any(is.na(diag(vcov.)))) {
            bad.se=T
        } else if (any(diag(vcov.)<0)) {
            bad.se=T
        } else {
            bad.se=F
        }            
        fit$bad.se=bad.se
        names(fit$coefficients)=substr(names(fit$coefficients), 1,1) # For now, it seems better to do this                
    }
    
    fit
} 
