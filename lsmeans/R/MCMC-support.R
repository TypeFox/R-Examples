# Support for MCMCglmm class and possibly more MCMC-based models

# Method to create a coda 'mcmc' object from a ref.grid
# (dots not supported, unfortunately)
as.mcmc.ref.grid = function(x, names = TRUE, ...) {
    object = x
    if (is.na(x@post.beta[1]))
        stop("No posterior sample -- can't make an 'mcmc' object")
    mat = x@post.beta %*% t(x@linfct)
    nm = setdiff(names(x@grid), c(".wgt.",".offset."))
    if (any(names)) {
        names = rep(names, length(nm))
        for (i in seq_along(nm))
            if(names[i]) x@grid[nm[i]] = paste(nm[i], x@grid[[nm[i]]])
    }
    dimnames(mat)[[2]] = do.call(paste, c(x@grid[, nm, drop = FALSE], sep=", "))
    coda::mcmc(mat, ...)
}

# Currently, data is required, as call is not stored
recover.data.MCMCglmm = function(object, data, ...) {    
    # if a multivariate response, stack the data with `trait` variable
    yvars = .all.vars(update(object$Fixed$formula, ". ~ 1"))
    if(length(yvars) > 1) {
#        for (v in yvars) data[[v]] = NULL
        dat = data
        for (i in seq_len(length(yvars) - 1))
            data = rbind(data, dat)
        data$trait = factor(rep(yvars, each = nrow(dat)))
    }
    attr(data, "call") = object$Fixed
    attr(data, "terms") = trms = delete.response(terms(object$Fixed$formula))
    attr(data, "predictors") = .all.vars(delete.response(trms))
    attr(data, "responses") = yvars
    data
}

lsm.basis.MCMCglmm = function(object, trms, xlev, grid, vcov., ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = NULL)
    Sol = as.matrix(object$Sol)[, seq_len(object$Fixed$nfl)] # toss out random effects if included
    bhat = apply(Sol, 2, mean)
    if (missing(vcov.))
        V = cov(Sol)
    else
        V = .my.vcov(object, vcov.)
    misc = list()
    list(X = X, bhat = bhat, nbasis = matrix(NA), V = V, 
         dffun = function(k, dfargs) NA, dfargs = list(), 
         misc = misc, post.beta = Sol)
}


### Support for MCMCpack , maybe others that produce mcmc objects
### Whether it works depends on:
###    1. if there is a "call" attribute with a formula or fixed member
###    2. if it's right, even then
### Alternatively, maybe providing formula and data will do the trick

recover.data.mcmc = function(object, formula, data, ...) {
    if (missing(formula)) {
        cl = attr(object, "call")
        if (is.null(cl$formula))
            cl$formula = cl$fixed
        if (is.null(cl$formula))
            return("No fixed-effects formula found")
        data = NULL
    }
    else {
        if (missing(formula) || missing(data))
            return("Requires both formula and data to proceed")
        cl = call("mcmc.proxy", formula = formula, data = quote(data))
    }
    trms = delete.response(terms(eval(cl$formula, parent.frame())))
    recover.data(cl, trms, NULL, data, ...)
}

lsm.basis.mcmc = function(object, trms, xlev, grid, vcov., ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = NULL)
    samp = as.matrix(object)[, seq_len(ncol(X)), drop = FALSE]
    bhat = apply(samp, 2, mean)
    if (missing(vcov.))
        V = cov(samp)
    else
        V = .my.vcov(object, vcov.)
    misc = list()
    list(X = X, bhat = bhat, nbasis = matrix(NA), V = V, 
         dffun = function(k, dfargs) NA, dfargs = list(), 
         misc = misc, post.beta = samp)
}


### support for CARBayes package - currently MUST supply data and have
### default contrasts matching what was used in fitting the mdoel
recover.data.carbayes = function(object, data, ...) {
    if(is.null(data)) # Try to recover data from parent frame
        data = model.frame(object$formula, data = parent.frame())
    cl = call("carbayes.proxy", formula = object$formula, data = quote(data))
    trms = delete.response(terms(eval(object$formula, parent.frame())))
    recover.data(cl, trms, NULL, data, ...)
}

lsm.basis.carbayes = function(object, trms, xlev, grid, ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = attr(object$X, "contrasts"))
    samp = as.matrix(object$samples$beta)
    bhat = apply(samp, 2, mean)
    V = cov(samp)
    misc = list()
    list(X = X, bhat = bhat, nbasis = matrix(NA), V = V, 
         dffun = function(k, dfargs) NA, dfargs = list(), 
         misc = misc, post.beta = samp)
}

