### Helper functions for lsmeans
### Here we have 'recover.data' and 'lsm.basis' methods
### For models that this package supports.
#----------------------------------------------------------
### Recover data methods will return a data.frame with 
### the original data, and at least these additional attrs:
#   attr(, "terms")      - terms component of object
#   attr(, "responses")  - names of response variable
#   attr(, "predictors") - names of predictors
#-----------------------------------------------------------
# generic version
recover.data = function(object, ...)
    UseMethod("recover.data")

#----------------------------------------------------------
### lsm.basis methods create a basis for the reference grid
#
# Required args:
#     object - the model object
#     trms   - terms component of object
#     xlev   - named list of factor levels (but not the coerced ones)
#     grid   - reference grid
# All methods must return a list with these elements:
#     X      - basis for linear fcns over reference grid
#     bhat   - regression coefficients for fixed effects (INCLUDING any NAs)
#     nbasis - matrix whose columns for a basis for non-estimable functions of beta; all.estble if none
#     V      - estimated covariance matrix of bhat
#     dffun  - function(k, dfargs) to find df for k'bhat having std error se
#     dfargs - additional arguments, if any, for dffun
#     misc   - Extra info ...
#              -- if extra levels need to be added (e.g. mlm, polr),
#                 put them in misc$ylevs
#              -- For transformations or link fcns, use misc$tran
#                 for name (see 'make.link'), and use misc$inv.lbl
#                 for label to use in 'summary' when tran is inverted
#                 (ref.grid looks at lhs of model for tran if none found)
# Note: if no df exists, set dffun = function(...) NA and dfargs = list()
#--------------------------------------------------------------
# generic version
lsm.basis = function(object, trms, xlev, grid, ...)
    UseMethod("lsm.basis")



#--------------------------------------------------------------
### DEFAULT METHODS (we hit these when a model is NOT supported)
# I'll have it return the message if we caught the error in this way
# Then caller can use try() to check for other types of errors,
# and just print this message otherwise 
recover.data.default = function(object, ...) {
    paste("Can't handle an object of class ", dQuote(class(object)[1]), "\n",
         paste(.show_supported(), collapse=""))
}

lsm.basis.default = function(object, trms, xlev, grid, ...) {
    stop("Can't handle an object of class", dQuote(class(object)[1]), "\n",
         .show_supported())
}

# Private fcn to get a list of supported objects
# does this by looking in namespace [ns] and methods [meth]
# then strips that off leaving extensions
.show_supported = function(ns = "lsmeans", meth = "lsm.basis") {
    "Use help(\"models\", package = \"lsmeans\") for information on supported models."
}


#--------------------------------------------------------------
### call' objects
# This recover.data method serves as the workhorse for the others
# For model objects, call this with the object's call and its terms component
# Late addition: if data is non-null, use it in place of recovered data
# Later addition: na.action arg req'd - vector of row indexes removed due to NAs
recover.data.call = function(object, trms, na.action, data = NULL, params = NULL, ...) {
    fcall = object # because I'm easily confused
    vars = setdiff(.all.vars(trms), params)
    if (length(vars) == 0)
        return("Model must have at least one predictor")
    tbl = data
    if (is.null(tbl)) {
        m = match(c("formula", "data", "subset", "weights"), names(fcall), 0L)
        fcall = fcall[c(1L, m)]
        fcall$drop.unused.levels = TRUE
        fcall[[1L]] = as.name("model.frame")
        fcall$xlev = NULL # we'll ignore xlev
        fcall$na.action = na.omit
        ### I once had a reason to put one var on the left, but don't remember why
        ### Now it messes up if there's a '$' in that term
        #         if (length(vars) > 1) 
        #             form = reformulate(vars[-1], response = vars[1])
        #         else
        form = reformulate(vars)
        fcall$formula = update(trms, form)
        env = environment(trms)
        if (is.null(env)) 
            env = parent.frame()
        tbl = eval(fcall, env, parent.frame())
        if (!is.null(na.action))
            tbl = tbl[-(na.action),  , drop=FALSE]
    }
    
    else
        fcall$data = tbl[complete.cases(data), , drop=FALSE]
    
    attr(tbl, "call") = object # the original call
    attr(tbl, "terms") = trms
    attr(tbl, "predictors") = setdiff(.all.vars(delete.response(trms)), params)
    attr(tbl, "responses") = setdiff(vars, union(attr(tbl, "predictors"), params))
    tbl
}


#--------------------------------------------------------------
### lm objects (and also aov, rlm, others that inherit) -- but NOT aovList
recover.data.lm = function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(terms(object)), object$na.action, ...)
}

lsm.basis.lm = function(object, trms, xlev, grid, ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    # coef() works right for lm but coef.aov tosses out NAs
    bhat = as.numeric(object$coefficients) 
    # stretches it out if multivariate - see mlm method
    V = .my.vcov(object, ...)
    
    if (sum(is.na(bhat)) > 0)
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    misc = list()
    if (inherits(object, "glm")) {
        misc = .std.link.labels(object$family, misc)
        dffun = function(k, dfargs) NA
        dfargs = list()
    }
    else {
        dfargs = list(df = object$df.residual)
        dffun = function(k, dfargs) dfargs$df
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
### mlm objects
# (recover.data.lm works just fine)

lsm.basis.mlm = function(object, trms, xlev, grid, ...) {
    class(object) = c("mlm", "lm") # avoids error in vcov for "maov" objects
    bas = lsm.basis.lm(object, trms, xlev, grid, ...)
    bhat = coef(object)
    k = ncol(bhat)
    bas$X = kronecker(diag(rep(1,k)), bas$X)
    bas$nbasis = kronecker(rep(1,k), bas$nbasis)
    ylevs = dimnames(bhat)[[2]]
    if (is.null(ylevs)) ylevs = seq_len(k)
    bas$misc$ylevs = list(rep.meas = ylevs)
    bas
}



#--------------------------------------------------------------
### merMod objects (lme4 package)
recover.data.merMod = function(object, ...) {
    if(!lme4::isLMM(object) && !lme4::isGLMM(object)) 
        return("Can't handle a nonlinear mixed model")
    fcall = object@call
    recover.data(fcall, delete.response(terms(object)), 
                 attr(object@frame, "na.action"), ...)
}

lsm.basis.merMod = function(object, trms, xlev, grid, vcov., ...) {
    if (missing(vcov.))
        V = as.matrix(vcov(object))
    else
        V = as.matrix(.my.vcov(object, vcov.))
    dfargs = misc = list()
    if (lme4::isLMM(object)) {
        pbdis = .lsm.is.true("disable.pbkrtest")
        Nlim = get.lsm.option("pbkrtest.limit")
        objN = lme4::getME(object, "N")
        toobig = objN > Nlim
        if (!pbdis && !toobig && requireNamespace("pbkrtest") && missing(vcov.)) {
            dfargs = list(unadjV = V, 
                adjV = pbkrtest::vcovAdj.lmerMod(object, 0))
            V = as.matrix(dfargs$adjV)
            tst = try(pbkrtest::Lb_ddf)
            if(class(tst) != "try-error")
                dffun = function(k, dfargs) pbkrtest::Lb_ddf (k, dfargs$unadjV, dfargs$adjV)
            else {
                dffun = function(k, dfargs) NA
                warning("To obtain d.f., install 'pbkrtest' version 0.4-1 or later")
            }
        }
        else {
            if(!pbdis && !("pbkrtest" %in% row.names(installed.packages())))
                message("Install package 'pbkrtest' to obtain bias corrections and degrees of freedom")
            else if(toobig)
                message("Note: Adjusted covariance and degrees-of-freedom calculations have been\n",
                        "disabled because the number of observations exceeds ", Nlim, ".\n",
                        "Standard errors and tests may be more biased than if they were adjusted.\n",
                        "To enable adjustments, set lsm.options(pbkrtest.limit = ", objN, ") or larger,\n",
                        "but be warned that this may result in large computation time and memory use.")
            dffun = function(k, dfargs) NA
        }
    }
    else if (lme4::isGLMM(object)) {
        dffun = function(k, dfargs) NA
        misc = .std.link.labels(family(object), misc)
    }
    else 
        stop("Can't handle a nonlinear mixed model")
    
    contrasts = attr(object@pp$X, "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = lme4::fixef(object)
    
    if (length(bhat) < ncol(X)) {
        # Newer versions of lmer can handle rank deficiency, but we need to do a couple of
        # backflips to put the pieces together right,
        # First, figure out which columns were retained
        kept = match(names(bhat), dimnames(X)[[2]])
        # Now re-do bhat with NAs in the right places
        bhat = NA * X[1, ]
        bhat[kept] = lme4::fixef(object)
        # we have to reconstruct the model matrix
        modmat = model.matrix(trms, object@frame, contrasts.arg=contrasts)
        nbasis = estimability::nonest.basis(modmat)
    }
    else
        nbasis=estimability::all.estble
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
### mer objects (from old lme4 version, now lme4.0)
### re-enabled; CRAN check now seems to work with multiple add'l repositories

recover.data.mer = function(object, ...) {
    if(!lme4.0::isLMM(object) && !lme4.0::isGLMM(object)) 
        return("Can't handle a nonlinear mixed model")
    fcall = object@call
    recover.data(fcall, delete.response(terms(object)), 
                 attr(object@frame, "na.action"), ...)
}

# Does NOT support pbkrtest capabilities. Uses asymptotic methods
lsm.basis.mer = function(object, trms, xlev, grid, ...) {
    V = as.matrix(.my.vcov(object, ...))
    dfargs = misc = list()
    if (lme4.0::isLMM(object)) {
        dffun = function(k, dfargs) NA        
    }
    else if (lme4.0::isGLMM(object)) {
        dffun = function(k, dfargs) NA
        # need to work harder as there is no 'family' method
        cfam = object@call$family
        if (is.name(cfam))
            fam = eval(cfam)()
        else
            fam = eval(cfam)
        misc = .std.link.labels(fam, misc)
    }
    else 
        stop("Can't handle a nonlinear mixed model")
    
    contrasts = attr(object@X, "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = lme4.0::fixef(object)
    nbasis=estimability::all.estble
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}




#--------------------------------------------------------------
### lme objects (nlme package)
recover.data.lme = function(object, data, ...) {
    fcall = object$call
    if (!is.null(fcall$weights))
        fcall$weights = nlme::varWeights(object$modelStruct)
    if(is.null(data)) # lme objects actually have the data, so use it!
        data = object$data
    recover.data(fcall, delete.response(terms(object)), object$na.action, data = data, ...)
}

lsm.basis.lme = function(object, trms, xlev, grid, sigmaAdjust = TRUE, ...) {
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = nlme::fixef(object)
    V = .my.vcov(object, ...)
    if (sigmaAdjust && object$method == "ML") 
        V = V * object$dims$N / (object$dims$N - nrow(V))
    misc = list()
    if (!is.null(object$family)) {
        misc = .std.link.labels(object$family, misc)
    }
    nbasis = estimability::all.estble
    # Replaced by containment method##dffun = function(...) NA
    dfx = object$fixDF$X
    if (names(bhat[1]) == "(Intercept)")
        dfx[1] = length(levels(object$groups[[1]])) - 1#min(dfx)   ### Correct apparent error in lme containment algorithm
    dffun = function(x, dfargs) {
        idx = which(abs(x) > 1e-4)
        ifelse(length(idx) > 0, min(dfargs$dfx[idx]), NA)
    }
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = list(dfx = dfx), misc = misc)
}



#--------------------------------------------------------------
### gls objects (nlme package)
recover.data.gls = function(object, ...) {
    fcall = object$call
    if (!is.null(fcall$weights))
        fcall$weights = nlme::varWeights(object$modelStruct)
    recover.data(fcall, delete.response(nlme::getCovariateFormula(object)), 
                 object$na.action, ...)
}

lsm.basis.gls = function(object, trms, xlev, grid, ...) {
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = coef(object)
    V = .my.vcov(object, ...)
    nbasis = estimability::all.estble
    dfargs = list(df = object$dims$N - object$dims$p)
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=list())
}



#--------------------------------------------------------------
### polr objects (MASS package)
recover.data.polr = recover.data.lm

lsm.basis.polr = function(object, trms, xlev, grid, 
                          mode = c("latent", "linear.predictor", "cum.prob", "prob", "mean.class"), 
                          rescale = c(0,1), ...) {
    mode = match.arg(mode)
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    # Strip out the intercept (borrowed code from predict.polr)
    xint = match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) 
        X = X[, -xint, drop = FALSE]
    bhat = c(coef(object), object$zeta)
    V = .my.vcov(object, ...)
    k = length(object$zeta)
    if (mode == "latent") {
        X = rescale[2] * cbind(X, matrix(- 1/k, nrow = nrow(X), ncol = k))
        bhat = c(coef(object), object$zeta - rescale[1] / rescale[2])
        misc = list(offset.mult = rescale[2])
    }
    else {
        j = matrix(1, nrow=k, ncol=1)
        J = matrix(1, nrow=nrow(X), ncol=1)
        X = cbind(kronecker(-j, X), kronecker(diag(1,k), J))
        link = object$method
        if (link == "logistic") link = "logit"
        misc = list(ylevs = list(cut = names(object$zeta)), 
                    tran = link, inv.lbl = "cumprob", offset.mult = -1)
        if (mode != "linear.predictor") {
            # just use the machinery we already have for the 'ordinal' package
            misc$mode = mode
            misc$postGridHook = ".clm.postGrid"
        }
    }
    misc$respName = as.character(terms(object))[2]
    nbasis = estimability::all.estble
    dffun = function(...) NA
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=list(), misc=misc)
}




#--------------------------------------------------------------
### survreg objects (survival package)
recover.data.survreg = function(object, ...) {
    fcall = object$call
    trms = delete.response(terms(object))
    # I'm gonna delete any terms involving strata(), cluster(), or frailty()
    mod.elts = dimnames(attr(trms, "factor"))[[2]]
    tmp = grep("strata\\(|cluster\\(|frailty\\(", mod.elts)
    if (length(tmp))
        trms = trms[-tmp]
    recover.data(fcall, trms, object$na.action, ...)
}

# Seems to work right in a little testing.
# However, it fails sometimes if I update the model 
# with a subset argument. Workaround: just fitting a new model
lsm.basis.survreg = function(object, trms, xlev, grid, ...) {
    # Much of this code is adapted from predict.survreg
    bhat = object$coefficients
    k = length(bhat)
    V = .my.vcov(object, ...)[seq_len(k), seq_len(k), drop=FALSE]
    # ??? not used... is.fixeds = (k == ncol(object$var))
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)    
    # X = model.matrix(object, m) # This is what predict.survreg does
    # But I have manipulated trms, so need to make sure things are consistent
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    nbasis = estimability::nonest.basis(model.matrix(object))
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df
    if (object$dist %in% c("exponential","weibull","loglogistic","loggaussian","lognormal")) 
        misc = list(tran = "log", inv.lbl = "response")
    else 
        misc = list()
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
###  coxph objects (survival package)
recover.data.coxph = function(object, ...) 
    recover.data.survreg(object, ...)

lsm.basis.coxph = function (object, trms, xlev, grid, ...) 
{
    object$dist = "doesn't matter"
    result = lsm.basis.survreg(object, trms, xlev, grid, ...)
    result$dfargs$df = NA
    result$X = result$X[, -1, drop = FALSE]
    result$X = result$X - rep(object$means, each = nrow(result$X))
    result$misc$tran = "log"
    result$misc$inv.lbl = "hazard"
    result
}

# Note: Very brief experimentation suggests coxph.penal also works.
# This is an extension of coxph


#--------------------------------------------------------------
###  coxme objects ####
### Greatly revised 6-15-15 (after version 2.18)
recover.data.coxme = function(object, ...) 
    recover.data.survreg(object, ...)

lsm.basis.coxme = function(object, trms, xlev, grid, ...) {
    bhat = fixef(object)
    k = length(bhat)
    V = .my.vcov(object, ...)[seq_len(k), seq_len(k), drop = FALSE]
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m)
    X = X[, -1, drop = FALSE] # remove the intercept
    # scale the linear predictor
    for (j in seq_along(X[1, ]))
        X[, j] = (X[, j] - object$means[j]) ### / object$scale[j]
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) NA
    misc = list(tran = "log", inv.lbl = "hazard")
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = list(), misc = misc)
}


###  special vcov prototype for cases where there are several vcov options
###  e.g., gee, geeglm, geese
.named.vcov = function(object, method, ...)
    UseMethod(".named.vcov")

# default has optional idx of same length as valid and if so, idx indicating 
#   which elt of valid to use if matched
# Ex: valid = c("mammal", "fish", "rat", "dog", "trout", "perch")
#     idx   = c(   1,        2,     1,     1,       2,       2)
#     -- so ultimately results can only be "mammal" or "fish"
# nonmatches revert to 1st elt.
.named.vcov.default = function(object, method, valid, idx = seq_along(valid), ...) {
    if (!is.character(method)) { # in case vcov. arg was matched by vcov.method {
        V = .my.vcov(object, method)
        method = "user-supplied"
    }
    else {
        i = pmatch(method, valid, 1)
        method = valid[idx[i]]
        V = object[[method]]
    }
    attr(V, "methMesg") = paste("Covariance estimate used:", method)
    V
}

# general-purpose lsm.basis function
.lsmb.geeGP = function(object, trms, xlev, grid, vcov.method, valid, idx = seq_along(valid), ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = coef(object)
    V = .named.vcov(object, vcov.method, valid, idx, ...)
    
    if (sum(is.na(bhat)) > 0)
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    
    misc = .std.link.labels(object$family, list())
    misc$initMesg = attr(V, "methMesg")
    dffun = function(k, dfargs) NA
    dfargs = list()
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}

#---------------------------------------------------------------
###  gee objects  ####


recover.data.gee = recover.data.lm

lsm.basis.gee = function(object, trms, xlev, grid, vcov.method = "robust.variance", ...)
    .lsmb.geeGP(object, trms, xlev, grid, vcov.method, 
                valid = c("robust.variance", "naive.variance"))

###  geepack objects  ####
recover.data.geeglm = recover.data.lm

lsm.basis.geeglm = function(object, trms, xlev, grid, vcov.method = "vbeta", ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = coef(object)
    V = .named.vcov(object$geese, vcov.method, 
                    valid = c("vbeta", "vbeta.naiv","vbeta.j1s","vbeta.fij","robust","naive"), 
                    idx = c(1,2,3,4,1,2))
    
    if (sum(is.na(bhat)) > 0)
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    
    misc = .std.link.labels(object$family, list())
    misc$initMesg = attr(V, "methMesg")
    dffun = function(k, dfargs) NA
    dfargs = list()
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}


recover.data.geese = function(object, ...) {
    fcall = object$call
    # what a pain - we need to reconstruct the terms component
    args = as.list(fcall[-1])
    na.action = object$na.action
    #trms = terms.formula(fcall$formula)
    if (!is.null(args$data)) {
        data = eval(args$data, parent.frame())
        trms = terms(model.frame(fcall$formula, data = data))
    } else {
        trms = terms(model.frame(fcall$formula))
    }
    recover.data(fcall, delete.response(trms), na.action, ...)
}

lsm.basis.geese = function(object, trms, xlev, grid, vcov.method = "vbeta", ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = object$beta
    V = .named.vcov(object, vcov.method, 
                    valid = c("vbeta", "vbeta.naiv","vbeta.j1s","vbeta.fij","robust","naive"), 
                    idx = c(1,2,3,4,1,2))

    # We don't have the qr component - I'm gonna punt for now
     if (sum(is.na(bhat)) > 0)
         warning("There are non-estimable functions, but estimability is NOT being checked")
#         nbasis = estimability::nonest.basis(object$qr)
#     else
        nbasis = estimability::all.estble
    
    misc = list()
    if (!is.null(fam <- object$call$family))
        misc = .std.link.labels(eval(fam)(), misc)
    misc$initMesg = attr(V, "methMesg")
    dffun = function(k, dfargs) NA
    dfargs = list()
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
### afex package - mixed objects ###
# just need to provide an 'lsmeans' method here, assuming Henrik adds the 'data' item

### These are deprecated as of afex 0.14 - now afex has its own lsmeans support

# recover.data.mixed = function(object, ...) {
#     recover.data.merMod(object$full.model, ...)
# }
# 
# lsm.basis.mixed = function(object, trms, xlev, grid, ...) {
#     lsm.basis.merMod(object$full.model, trms, xlev, grid, ...)
# }


#--------------------------------------------------------------
### glmmADMB package

recover.data.glmmadmb = recover.data.lm

lsm.basis.glmmadmb = function (object, trms, xlev, grid, ...) 
{
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = glmmADMB::fixef(object)
    V = .my.vcov(object, ...)
    misc = list()
    if (!is.null(object$family)) {
        fam = object$family
        misc$tran = object$link
        misc$inv.lbl = "response"
        if (!is.na(pmatch(fam,"binomial"))) 
            misc$inv.lbl = "prob"
        else if (!is.na(pmatch(fam,"poisson"))) 
            misc$inv.lbl = "rate"
    }
    nbasis = estimability::all.estble
    dffun = function(...) NA
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = list(), misc = misc)
}


# --------------------------------------------------------------
### Explicit non-support for 'gam' objects (runs, but results are wrong)

lsm.basis.gam = function(object, trms, xlev, grid, ...) {
    stop("Can't handle an object of class", dQuote(class(object)[1]), "\n",
         .show_supported())
}


#--------------------------------------------------------------
### mgcv package --------------------------------

# gam - OK - inherits from glm

# gamm
# recover.data.gamm = function(object, ...) {
#     fcall = object$lme$call
#     recover.data(fcall, delete.response(terms(object$lme)), object$lme$na.action, ...)
# }
# 
# lsm.basis.gamm = function (object, trms, xlev, grid, sigmaAdjust = TRUE, ...) {
#     lsm.basis(object$lme, trms, xlev, grid, sigmaAdjust, ...)
#     # Doesn't work because needs the matrices in object$lme$data
# }






### ----- Auxiliary routines -------------------------
# Provide for vcov. argument in ref.grid call, which could be a function or a matrix

.my.vcov = function(object, vcov. = stats::vcov, ...) {
    if (is.function(vcov.))
        vcov. = vcov.(object)
    else if (!is.matrix(vcov.))
        stop("vcov. must be a function or a square matrix")
    vcov.
}

# Call this to do the standard stuff with link labels
# Returns a modified misc
.std.link.labels = function(fam, misc) {
    if (is.null(fam))
        return(misc)
    misc$tran = fam$link
    misc$inv.lbl = "response"
    if (length(grep("binomial", fam$family)) == 1)
        misc$inv.lbl = "prob"
    else if (length(grep("poisson", fam$family)) == 1)
        misc$inv.lbl = "rate"
    misc
}

## Alternative to all.vars, but keeps vars like foo$x and foo[[1]] as-is
##   Passes ... to all.vars
.all.vars = function(expr, retain = c("\\$", "\\[\\[", "\\]\\]"), ...) {
    if (!inherits(expr, "formula")) {
        expr = try(eval(expr), silent = TRUE)
        if(inherits(expr, "try-error")) {
            return(character(0))
        }
    }
    repl = paste("_Av", seq_along(retain), "_", sep = "")
    for (i in seq_along(retain))
        expr = gsub(retain[i], repl[i], expr)
    subs = switch(length(expr), 1, c(1,2), c(2,1,3))
    vars = all.vars(as.formula(paste(expr[subs], collapse = "")), ...)
    retain = gsub("\\\\", "", retain)
    for (i in seq_along(retain))
        vars = gsub(repl[i], retain[i], vars)
    vars
}


### Not-so-damn-smart replacement of diag() that will 
### not be so quick to assume I want an identity matrix
### returns matrix(x) when x is a scalar
.diag = function(x, nrow, ncol) {
    if(is.matrix(x))
        diag(x)
    else if((length(x) == 1) && missing(nrow) && missing(ncol)) 
        matrix(x)
    else 
        diag(x, nrow, ncol)
}


