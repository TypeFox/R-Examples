# effect.mer and effect.lme built from effect.lm by S. Weisberg 29 June 2011
# last modified 2012-03-08 to require() lme4 or nlme. J. Fox
# 2012-10-05 effect.lme didn't work with 'weights', now corrected.  S. Weisberg
# 2013-03-05: introduced merMod methods for development version of lme4. J. Fox
# 2013-04-06: added support for lme4.0, J. Fox
# 2013-07-30: added 'data' argument to lme.to.glm and mer.to.glm to allow
#   calling effect from within a subroutine.
# 2013-09-25:  removed the 'data' argument as it make the functions fail with
#   logs, splines and polynomials
# 2014-09-24: added option for KR cov matrix to mer.to.glm(). J. Fox
# 2014-12-07: don't assume that pbkrtest is installed. J. Fox
# 2014-12-20: mer.to.glm failed for negative.binomial() because the link has an argument
#   that was handled incorrectly by the family.glmResp function.  This function is no longer
#   used by mer.to.glm.  The same error will recur in any link with an argument.
# 2015-06-10: requireNamespace("pbkrtest") rather than require("pbkrtest)
# 2015-07-02: fixed bug when the name of the data frame was the name of a function (e.g., sort, or lm)
# 2015-12-13: make it work with pbkrtest 0.4-3. J. Fox
# 2016-01-07: modified 'fixmod' to allow "||" in variance formulae
# 2016-01-19: Fixed bug in glm.to.mer when 'poly' is used in a model.

# the function lm.wfit fit gets the hessian wrong for mer's.  Get the variance
# from the vcov method applied to the mer object.


fixmod <- function (term) 
{
  if (!("|" %in% all.names(term)) && !("||" %in% all.names(term)))
    return(term)
  if ((is.call(term) && term[[1]] == as.name("|")) ||
      (is.call(term) && term[[1]] == as.name("||")))
    return(NULL)
    if (length(term) == 2) {
        nb <- fixmod(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- fixmod(term[[2]])
    nb3 <- fixmod(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

# lme.to.glm evaluates a 'glm' model that is as similar to a given 'lme'
# model, in the same pattern as mer.to.glm.  This could be speeded up
# slightly by using 'lm' rather than 'glm' but I use 'glm' to parallel
# mer.to.glm more closely.  The differences are:  (1) match fewer args
# in the call; (2) different def of mod2$coefficients; no other 
# changes

# The argument 'data' to lme.to.glm and to mer.to.glm copies the data
# from the object into the local environment and makes it visible when 'effect'
# is called from within another function.

lme.to.glm <- function(mod) {
    cl <- mod$call
    cl$formula <- cl$fixed
    m <- match(c("formula", "data", "subset", 
        "na.action",  "contrasts"), names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- as.name("glm")
    mod2 <- eval(cl)
    pw <- attr(mod$modelStruct$varStruct, "weights")
    if(!is.null(pw)) mod2$prior.weights <- pw
    mod2$coefficients <- mod$coefficients$fixed
    mod2$vcov <- as.matrix(vcov(mod))
    mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
    mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
    mod2$weights <- as.vector(with(mod2,
        prior.weights * (family$mu.eta(linear.predictors)^2 /
                family$variance(fitted.values))))
    mod2$residuals <- with(mod2,
        prior.weights * (y - fitted.values)/weights )
    class(mod2) <- c("fakeglm", class(mod2))
    mod2
}

# mer.to.glm evaluates a 'glm' model that is as similar to a given 'mer'
# model as follows.  It is of class c("fakeglm", "glm", "lm")
# several items are added to the created objects. Do not export

mer.to.glm <- function(mod, KR=FALSE) {
    if (KR && !requireNamespace("pbkrtest", quietly=TRUE)){
        KR <- FALSE
        warning("pbkrtest is not available, KR set to FALSE")
    }
    # object$family$family doesn't work correctly with the negative binomial family because of the
    # argument in the family function, so the old line
    #   family <- family(mod)
    # returns an error message for these models.  The following kluge fixes this.
    # If this bug is fixed in lme4, this code may break because it expects resp$family$family
    # to return "Link Name(arg)" with ONE argument, and so spaces between Name and "(arg)"
    family1 <- function(object, ...) {UseMethod("family1", object@resp)}
    family1.lmResp <- function(object, ...) family(object, ...)
    family1.glmResp <- function(object, ...){
        famname <- object@resp$family$family
        open.paren <- regexpr("\\(", famname)
        if(open.paren==-1) {
            name <- famname
            arg <- list()
        } else {
            name <- sub(" ", ".", tolower(substr(famname, 1, -1 + open.paren)))
            arg <- list(as.numeric(gsub("\\)", "", substr(famname, 1 + open.paren, 100))))
        }
        if(is.null(object@resp$family$initialize)) do.call(name, arg) else object@resp$family
    }
    family <- family1(mod)
    # end
    link <- family$link
    family <- family$family
    cl <- mod@call
    if(cl[[1]] =="nlmer") stop("effects package does not support 'nlmer' objects")
    m <- match(c("formula", "family", "data", "weights", "subset", 
        "na.action", "start", "offset",  
        "model", "contrasts"), names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- as.name("glm")
    cl$formula <- fixmod(as.formula(cl$formula))
#    cl$data <- mod@frame # caused bug with a 'poly' in the formula
    mod2 <- eval(cl)
    mod2$coefficients <- lme4::fixef(mod) #mod@fixef
    mod2$vcov <- if (family == "gaussian" && link == "identity" && KR) as.matrix(pbkrtest::vcovAdj(mod)) else as.matrix(vcov(mod))
    mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
    mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
    mod2$weights <- as.vector(with(mod2,
        prior.weights * (family$mu.eta(linear.predictors)^2 /
                family$variance(fitted.values))))
    mod2$residuals <- with(mod2,
        prior.weights * (y - fitted.values)/weights )
    class(mod2) <- c("fakeglm", class(mod2))
    mod2
}



#method for 'fakeglm' objects. Do not export   
vcov.fakeglm <- function(object, ...) object$vcov

#The next six functions should be exported as S3 methods

effect.mer <- function(term, mod, vcov.=vcov, KR=FALSE, ...) {
    result <- effect(term, mer.to.glm(mod, KR=KR), vcov., ...)
    result$formula <- as.formula(formula(mod))
    result
}

effect.merMod <- function(term, mod, vcov.=vcov, KR=FALSE, ...){
    effect.mer(term, mod, vcov.=vcov, KR=KR, ...)
}

effect.lme <- function(term, mod, ...) {
    mod1 <- lme.to.glm(mod)
    result <- effect(term, mod1)
    result$formula <- as.formula(formula(mod))
    result
}

allEffects.mer <- function(mod, KR=FALSE,...){
    allEffects(mer.to.glm(mod,KR=KR), ...)
}

allEffects.merMod <- function(mod, KR=FALSE,...){
    allEffects(mer.to.glm(mod,KR=KR), ...)
}

allEffects.lme <- function(mod, ...){
    allEffects(lme.to.glm(mod), ...)
}



