
# $Id: helpers.R 412 2015-02-03 17:19:12Z thothorn $

### model.matrix.coxph doesn't return contrasts etc.
model.matrix.coxph <- function(object, ...) {
    mm <- model.matrix(delete.response(terms(object)),
                       data = model.frame(object))
    at <- attributes(mm)
    mm <- mm[,-1]
    at$dim[2] <- at$dim[2] - 1
    at$dimnames[[2]] <- at$dimnames[[2]][-1]
    at$assign <- at$assign[-1]
    attributes(mm) <- at
    mm
}

model.matrix.coxph.penal <- function(object, ...) {

    class(object) <- "coxph"
    mm <- model.matrix(object)
    at <- attributes(mm)
    indx <- grep("frailty", colnames(mm))
    ret <- mm[ , -indx, drop = FALSE]
    attr(ret, "assign") <- at$assign[-indx]
    attr(ret, "contrasts") <- at$contrasts
    ret
}

model.frame.coxph.penal <- function(object, ...) {

    tm <- terms(object)
    class(object) <- "coxph"
    mf <- model.frame(object)
    ret <- cbind(mf[[1]], model.frame(delete.response(tm), data = mf))
    colnames(ret)[1] <- colnames(mf)[1]
    ret
}

terms.coxph.penal <- function(object, ...) {

    class(object) <- "coxph"           
    tm <- terms(object)
    ctm <- as.character(tm)
    x <- strsplit(ctm[3], "+", fixed = TRUE)[[1]]
    x <- x[-grep("frailty", x)]
    fm <- paste(ctm[2], "~", paste(x, collapse = "+"))
    terms(as.formula(fm))
}

coxph.penalcoef <- function(object, ...) {

    mm <- model.matrix(object)
    class(object) <- "coxph"
    cf <- coef(object)
    cf[1:ncol(mm)]
}

coxph.penalvcov <- function(object, ...) {
    
    mm <- model.matrix(object)
    class(object) <- "coxph"
    vc <- vcov(object)        
    vc[1:ncol(mm), 1:ncol(mm), drop = FALSE]
}


model.matrix.survreg <- function(object, ...) {
   model.matrix(delete.response(terms(object)),
                       data = model.frame(object))
}

### coxme objects
model.matrix.coxme <- model.matrix.coxph



model.matrix.aovlist <- function(object, ...)
    stop(sQuote("glht"), " does not support objects of class ", 
         sQuote("aovlist"))

model.matrix.lme <- function(object, ...)
    model.matrix(terms(object), data = model.frame(object), ...)

model.frame.lme <- function(object, ...) {
    ret <- object$data
    if (is.null(ret)) stop("object does not contain any data")
    ret
}

### extract coefficients, covariance matrix and 
### degrees of freedom (if available) from `model'
modelparm <- function(model, coef., vcov., df, ...) 
    UseMethod("modelparm")

modelparm.default <- function(model, coef. = coef, vcov. = vcov, 
                              df = NULL, ...) 
{

    ### allow specification of coef and vcov directly
    if (!is.function(coef.)) {
        beta <- coef.
        coef. <- function(model) return(beta)
    }
    if (!is.function(vcov.)) {
        sigma <- vcov.
        vcov. <- function(model) return(sigma)
    }

    ### extract coefficients and their covariance matrix
    beta <- try(coef.(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    sigma <- try(vcov.(model))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("model"), " found!")       
    sigma <- as.matrix(sigma)

    if (any(length(beta) != dim(sigma))) 
        stop("dimensions of coefficients and covariance matrix don't match")

    ### determine degrees of freedom
    if (is.null(df)) {
        df <- 0
        ### check if a linear model was supplied
        if (class(model)[1] %in% c("aov", "lm")) {
            class(model) <- "lm"
            df <- summary(model)$df[2]
        }
        if (inherits(model, "parm"))
            df <- model$df
    } else {
        if (df < 0) stop(sQuote("df"), " is not positive")
    }

    ### try to identify non-estimable coefficients
    ### coef.aov removes NAs, thus touch coefficients 
    ### directly
    ocoef <- coef.(model)
    if (inherits(model, "aov")) ocoef <- model$coefficients
    estimable <- rep(TRUE, length(ocoef))
    if (any(is.na(ocoef))) {
        estimable[is.na(ocoef)] <- FALSE
        beta <- ocoef[estimable]
    }

    ### just in case...
    if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
        stop("could not extract coefficients and covariance matrix from ", 
             sQuote("model"))

    RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
    class(RET) <- "modelparm"
    RET
}

### mixed effects models (package `lme4')
modelparm.mer <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### mixed effects models (package `lme4Eigen')
modelparm.merMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### package `nlme'
modelparm.lme <- function(model, coef. = nlme::fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### survreg models (package `survival')
vcovsurvreg <- function(object, ...) {
    sigma <- vcov(object)
    p <- length(coef(object))
    return(sigma[1:p, 1:p])
}

modelparm.survreg <- function(model, coef. = coef, vcov. = vcovsurvreg, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

modelparm.aovlist <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    stop(sQuote("glht"), " does not support objects of class ", sQuote("aovlist"))

modelparm.coxme <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

modelparm.coxph.penal <- function(model, coef. = coxph.penalcoef, 
                                  vcov. = coxph.penalvcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

model.matrix.polr <- model.matrix.coxph

polrvcov <- function(object) {
   cf <- coef(object)
   vcov <- vcov(object)
   vcov[names(cf), names(cf)]
}

modelparm.polr <- function(model, coef. = coef, vcov. = polrvcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### modified from package MASS  
MPinv <- function (X, tol = sqrt(.Machine$double.eps))
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) 
        RET <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))   
    else if (!any(Positive))
        RET <- array(0, dim(X)[2:1])
    else RET <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
        t(Xsvd$u[, Positive, drop = FALSE]))
    return(list(MPinv = RET, rank = sum(Positive)))
}

### meaningless ...
chkdots <- function(...) {

    lst <- list(...)
    if (length(lst) > 0) {
        warning("Argument(s) ", sQuote(names(lst)), " passed to ", sQuote("..."), 
                " are ignored", call. = TRUE)
    }
}
