pmpmodel <-
function (bmao, model = numeric(0), exact = TRUE) 
{
    if (!is.bma(bmao)) 
        stop("bmao needs to be a bma object")
    if (!is.vector(model)) 
        stop("model needs to be vector denoting a single model.")
    K = bmao$info$K
    was.enum = (bmao$arguments$mcmc == "enum")
    emptyindex = logical(K)
    modelhex = ""
    if (length(model) == 0L) 
        model = numeric(0)
    if ((is.character(model)) && (all(model %in% bmao$reg.names))) {
        mix = match(model, bmao$reg.names)
        if (any(is.na(mix))) 
            stop("Provided variable names do not conform to bma object")
        emptyindex[mix] = TRUE
        model = emptyindex
    }
    else if ((length(model) == 1L) && all(is.character(model))) {
        modelhex = model
        model = as.logical(hex2bin(model))
        if (length(model) > K) 
            model = model[-(1:(length(model) - K))]
    }
    else if (is.logical(model) || ((length(model) == K) && (is.numeric(model) && 
        max(model) < 2))) {
        if (length(model) > K) 
            model = model[-(1:(length(model) - K))]
        model = as.logical(model)
    }
    else if (is.numeric(model)) {
        emptyindex[model] = TRUE
        model = emptyindex
    }
    else stop("model needs to be an integer, logical or character model index representation (hexcode or variable names)")
    if (any(is.na(model))) 
        stop("Provided model index does not seem to exist.")
    if (modelhex == "") 
        modelhex = bin2hex(model)
    fixed.pos = bmao$mprior.info$fixed.pos
    if (is.null(fixed.pos)) 
        fixed.pos = numeric(0)
    if (any(model[fixed.pos] != TRUE)) 
        stop("Such a model was excluded by bmao's argument fixed.reg")
    bools = bmao$topmod$bool()
    liks = bmao$topmod$lik()
    ncounts = bmao$topmod$ncount()
    cumsumweights = bmao$info$cumsumweights
    yty = as.vector(crossprod(bmao$X.data[, 1, drop = TRUE] - 
        mean(bmao$X.data[, 1, drop = TRUE])))
    log.null.lik = bmao$gprior.info$lprobcalc$just.loglik(ymy = yty, 
        k = 0)
    mix = match(modelhex, bools)
    if ((!exact) && (!was.enum)) {
        if (!is.na(mix)) {
            return(ncounts[[mix]]/cumsumweights)
        }
        else if (!any(model) || all(model)) {
            return(bmao$info$k.vec[sum(model) + 1]/cumsumweights)
        }
        else {
            stop("Model MCMC-based PMP cannot be found. Try exact=TRUE .")
        }
    }
    if (!is.na(mix)) {
        loglik = liks[mix]
    }
    else if (was.enum && (!any(model) || all(model))) {
        loglik = log(bmao$info$k.vec[sum(model) + 1]) + log.null.lik
    }
    else {
        if (!was.enum && (length(liks) == 0L)) 
            stop("bmao needs to contain more than 0 top models to provide an estimate for your model index.")
        if (sum(model) == 0L) {
            loglik = log.null.lik + bmao$mprior.info$pmp(ki = 0, 
                mdraw = rep(0, K), ymy = yty)
        }
        else {
            zz = zlm(bmao$X.data[, c(TRUE, model), drop = FALSE], 
                g = bmao$gprior.info)
            loglik = zz$marg.lik + bmao$mprior.info$pmp(ki = sum(model), 
                mdraw = as.numeric(model), ymy = zz$olsres$ymy)
        }
    }
    if (was.enum) {
        return(exp(loglik - log.null.lik)/cumsumweights)
    }
    pmp_withintop = exp(loglik - log.null.lik)/sum(exp(liks - 
        log.null.lik))
    return(pmp_withintop * sum(ncounts)/cumsumweights)
}
