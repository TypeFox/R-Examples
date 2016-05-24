# experimental support for nls, nlme objects

recover.data.nls = function(object, ...) {
    fcall = object$call
    trms = terms(reformulate(names(object$dataClasses)))
    recover.data(fcall, trms, object$na.action, ...)
}

lsm.basis.nls = function(object, trms, xlev, grid, ...) {
    Vbeta = .my.vcov(object, ...)
    env = object$m$getEnv()
    for (nm in names(grid)) env[[nm]] = grid[[nm]]
    pars = object$m$getAllPars()
    DD = deriv(object$m$formula(), names(pars))
    ests = eval(DD, env)
    bhat = as.numeric(ests)
    grad = attr(ests, "gradient")
    V = grad %*% Vbeta %*% t(grad)
    X = diag(1, nrow(grid))
    list(X=X, bhat=bhat, nbasis=all.estble, V=V, 
         dffun=function(k, dfargs) NA, dfargs=list(), 
         misc=list())
}
    

### For nlme objects, we can do stuff with the fixed part of the model
### Additional REQUIRED argument is 'param' - parameter name to explore
recover.data.nlme = function(object, param, ...) {
    if(missing(param))
        return("'param' argument is required for nlme objects")
    fcall = object$call
    if (!is.null(fcall$weights))
        fcall$weights = nlme::varWeights(object$modelStruct)
    fixed = fcall$fixed
    if (is.call(fixed))
        fixed = eval(fixed, envir = parent.frame())
    if(!is.list(fixed))
        fixed = list(fixed)
    form = NULL
    for (x in fixed)
        if (param %in% all.names(x)) form = x
    if (is.null(form))
        return(paste("Can't find '", param, "' among the fixed parameters", sep = ""))
    fcall$weights = NULL
    #trms = delete.response(terms(update(terms(object), form)))
    trms = delete.response(terms(form))
    if (length(.all.vars(trms)) == 0)
        return(paste("No predictors for '", param, "' in fixed model", sep = ""))
    recover.data(fcall, trms, object$na.action, ...)
}

lsm.basis.nlme = function(object, trms, xlev, grid, param, ...) {
    idx = object$map$fmap[[param]]
    V = object$varFix[idx, idx, drop = FALSE]
    bhat = object$coefficients$fixed[idx]
    contr = attr(object$plist[[param]]$fixed, "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contr)
    dfx = object$fixDF$X[idx]
    dfx[1] = min(dfx) # I'm assuming 1st one is intercept
    dffun = function(k, dfargs) { # containment df
        idx = which(abs(k) > 1e-6)
        ifelse(length(idx) > 0, min(dfargs$dfx[idx]), NA)
    }
    list(X = X, bhat = bhat, nbasis = estimability::all.estble, 
         V = V, dffun = dffun, dfargs = list(dfx = dfx), 
         misc = list(estName = param))
}


