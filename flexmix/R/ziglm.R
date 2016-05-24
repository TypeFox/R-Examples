#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: ziglm.R 4834 2012-08-02 10:17:09Z gruen $
#

setClass("FLXMRziglm", contains = "FLXMRglm")

FLXMRziglm <- function(formula = . ~ .,
                       family = c("binomial", "poisson"), ...) {
  family <- match.arg(family)
  new("FLXMRziglm", FLXMRglm(formula, family, ...),
      name = paste("FLXMRziglm", family, sep=":"))
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRziglm"),
          function(model, data, formula, lhs=TRUE, ...) {
  model <- callNextMethod(model, data, formula, lhs)
  if (attr(terms(model@fullformula), "intercept") == 0)
    stop("please include an intercept")
  model
})

setMethod("FLXremoveComponent", signature(model = "FLXMRziglm"),
          function(model, nok, ...) 
          if (1 %in% nok) as(model, "FLXMRglm") else model)

setMethod("FLXmstep", signature(model = "FLXMRziglm"),
          function(model, weights, components, ...) {
  coef <- c(-Inf, rep(0, ncol(model@x)-1))
  names(coef) <- colnames(model@x)
  comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                 family = model@family), eval(model@defineComponent))
  c(list(comp.1),
    FLXmstep(as(model, "FLXMRglm"), weights[, -1, drop=FALSE], components[-1]))
})

setMethod("FLXgetDesign", signature(object = "FLXMRziglm"),
function(object, components)
  rbind(0, FLXgetDesign(as(object, "FLXMRglm"), components[-1])))

setMethod("FLXreplaceParameters", signature(object="FLXMRziglm"),
function(object, components, parms)
          c(components[[1]], FLXreplaceParameters(as(object,
            "FLXMRglm"), components[-1], parms)))

setMethod("FLXgradlogLikfun", signature(object="FLXMRziglm"),
function(object, components, weights, ...)
          FLXgradlogLikfun(as(object, "FLXMRglm"),
                           components[-1], weights[,-1,drop=FALSE]))

setMethod("refit_optim", signature(object = "FLXMRziglm"),
function(object, components, ...) {
  x <- refit_optim(as(object, "FLXMRglm"), components[-1], ...)
  names(x) <- paste("Comp", 1 + seq_along(x), sep = ".")
  x
})
