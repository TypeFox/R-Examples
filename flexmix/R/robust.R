#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: robust.R 4834 2012-08-02 10:17:09Z gruen $
#
###*********************************************************

setClass("FLXMRrobglm",
         representation(bgw="logical"),
         prototype(bgw=FALSE),
         contains = "FLXMRglm")


FLXMRrobglm <- function(formula = . ~ .,
                      family=c("gaussian", "poisson"),
                      bgw=FALSE, ...)
{
    family <- match.arg(family)
    new("FLXMRrobglm", FLXMRglm(formula, family, ...),
        name = paste("FLXMRrobglm", family, sep=":"),
        bgw = bgw)
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRrobglm"),
function(model, data, formula, lhs=TRUE, ...)
{
    model <- callNextMethod(model, data, formula, lhs)
    if (attr(terms(model@fullformula), "intercept")==0)
        stop("please include an intercept")
    new("FLXMRrobglm", model)
})

setMethod("FLXremoveComponent", signature(model = "FLXMRrobglm"),
          function(model, nok, ...)
{
    if (1 %in% nok) model <- as(model, "FLXMRglm")
    model
})

setMethod("FLXmstep", signature(model = "FLXMRrobglm"),
          function(model, weights, ...)
{
    if(model@bgw){
        w <- weights[,1]
    }
    else{
        w <- rep(1, nrow(weights))
    }

    if(model@family=="gaussian")
    {
        cwt <- cov.wt(model@y, w)
        coef <- c(cwt$center, rep(0, ncol(model@x)-1))
    
        names(coef) <- colnames(model@x)
        comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                            sigma=sqrt(cwt$cov),
                            family = model@family),
                       eval(model@defineComponent))
    }
    else if(model@family=="poisson")
    {
        cwt <- cov.wt(model@y, w)
        coef <- c(log(3*cwt$center), rep(0, ncol(model@x)-1))
        names(coef) <- colnames(model@x)
        
        comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                            family = model@family),
                       eval(model@defineComponent))
    }
    else{
        stop("Other families not implemented yet!")
    }
    
    c(list(comp.1), FLXmstep(as(model, "FLXMRglm"),
                             weights[, -1, drop=FALSE]))
})

