# if this file is sourced during development, also need to source ssfct.R

# it seems that drm expects weights to be a global variable, not formula, not local variable, although it evaluates weights in parent.frame(), so I don't really understand
# r cmd check does not allow global variables 
if(getRversion() >= "2.15.1")  utils::globalVariables(c("drm.weights"))

# robust="mean"; fit.4pl=FALSE; weighting=FALSE # default
drm.fit=function (formula, data, robust="mean", fit.4pl=FALSE, w=NULL, gof.threshold=0.2, verbose=FALSE, bcVal = NULL, bcAdd = 0) {
    
    force.fit=TRUE
    # try three self starting functions
    gofs=rep(NA, 3)
    fits=list()
    outcome.coln=all.vars(formula)[1]
    predictor.coln=all.vars(formula)[2]
    
    if (is.null(w)) {
        # using assign function with inherits, instead of <<-, avoids a note during r cmd check. 
        # using assign function with .GlobalEnv still gets a not by r cmd check
        # if rcmdcheck notes the following, it will be reall hard to find a way around
        assign("drm.weights", rep(1,nrow(data)), inherits=TRUE)
        #drm.weights <<- rep(1,nrow(data))
    } else {
        stopifnot(length(w)==nrow(data))
        assign("drm.weights", w, inherits=TRUE)
        #drm.weights <<- w
        stopifnot(length(drm.weights)==nrow(data))
    }
        
    control=drmc(maxIt=500, method="BFGS", relTol=1e-7, trace=FALSE)

    bad.se=TRUE 
    if (!fit.4pl){
        # try default ss.fct first so that the results will match drm call in most cases
        fit1= try(drm(formula=formula, data=data, robust=robust, fct = LL.5(), control=control, weights=drm.weights, bcVal=bcVal, bcAdd=bcAdd), silent=verbose>=2)    
        
        if (!inherits(fit1, "try-error")) {
            vcov.=vcov(fit1)
            if(is.null(vcov.)) {
                bad.se = TRUE 
            } else if(any(is.na(diag(vcov.)))) {
                bad.se=TRUE
            } else if (any(diag(vcov.)<0)) {
                bad.se=TRUE
            } else {
                bad.se=FALSE
            }            
            gof1 = 1-neill.test(fit1,display=FALSE) # 1-p val, when display is FALSE, F-value is not returned
        } else {
            gof1=Inf
        }
        gofs[1]=gof1
        fits[[1]]=fit1
        
        if (1-gof1<gof.threshold) {
            
            if (verbose) cat("drm.fit: default ssfct fails, try additional starting functions.  ")
            fit2= try(drm(formula=formula, data=data, robust=robust, fct = LL.5(ssfct=ssfct.drc.1.5.2), control=control, weights=drm.weights, bcVal=bcVal, bcAdd=bcAdd), silent=verbose>=2)
            
            if (!inherits(fit2, "try-error")) {
                vcov.=vcov(fit2)
                if(is.null(vcov.)) {
                    bad.se = TRUE 
                } else if(all(is.na(diag(vcov.)))) {
                    bad.se=TRUE
                } else if (any(diag(vcov.)<0)) {
                    bad.se=TRUE
                } else {
                    bad.se=FALSE
                }            
                gof2 = 1-neill.test(fit2,display=FALSE)
            } else gof2=Inf
            gofs[2]=gof2
            fits[[2]]=fit2
            
            fit3 = try(drm(formula=formula, data=data, robust=robust, fct = LL.5(ssfct=ss.fct.via.LL4), control=control, weights=drm.weights, bcVal=bcVal, bcAdd=bcAdd), silent=verbose>=2)
            
            if (!inherits(fit3, "try-error")) {
                vcov.=vcov(fit3)
                if(is.null(vcov.)) {
                    bad.se = TRUE 
                } else if(all(is.na(diag(vcov.)))) {
                    bad.se=TRUE
                } else if (any(diag(vcov.)<0)) {
                    bad.se=TRUE
                } else {
                    bad.se=FALSE
                }            
                gof3 = 1-neill.test(fit3,display=FALSE)
            } else gof3=Inf
            gofs[3]=gof3
            fits[[3]]=fit3
                        
        } 
        
        fit=fits[[which.min(gofs)]]
        
    } else {
        fit= try(drm(formula=formula, data=data, robust=robust, fct = LL.4(), control=control, weights=drm.weights, bcVal=bcVal, bcAdd=bcAdd), silent=verbose>=2)  
        if (!inherits(fit, "try-error")) {
            vcov.=vcov(fit)
            if(is.null(vcov.)) {
                bad.se = TRUE 
            } else if(all(is.na(diag(vcov.)))) {
                bad.se=TRUE
            } else if (any(diag(vcov.)<0)) {
                bad.se=TRUE
            } else {
                bad.se=FALSE
            }            
            if (bad.se) gof3=Inf
            else gof3 = 1-neill.test(fit,display=FALSE)
        } else gof3=Inf
        gofs[1]=gof3
        
    }
    
    if (verbose) myprint(gofs)
    if (class(fit)[1]=="try-error") {
        if (verbose) cat("return null fit\n")
        return (NULL)
    } else {
        fit$bad.se=bad.se
        names(fit$coefficients)=substr(names(fit$coefficients), 1,1) # For now, it seems better to do this
        return (fit)
    }
    
}



getVarComponent.drc=function (object,...) {
    summary(object)$resVar
}
