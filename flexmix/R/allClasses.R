#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: allClasses.R 5016 2015-01-17 11:12:41Z gruen $
#

setClass("FLXcontrol",
         representation(iter.max="numeric",
                        minprior="numeric",
                        tolerance="numeric",
                        verbose="numeric",
                        classify="character",
                        nrep="numeric"),
         prototype(iter.max=200,
                   minprior=0.05,
                   tolerance=10e-7,
                   verbose=0,
                   classify="auto",
                   nrep=1),
         validity=function(object) {
           (object@iter.max > 0)
         })

setAs("list", "FLXcontrol",
function(from, to){
    z = list2object(from, to)
    z@classify = match.arg(z@classify,
                           c("auto", "weighted", "hard", "random", "SEM", "CEM"))
    z
})

setAs("NULL", "FLXcontrol",
function(from, to){
    new(to)
})


###**********************************************************

setClass("FLXM",
         representation(fit="function",
                        defineComponent="expression",
                        weighted="logical",
                        name="character",
                        formula="formula",
                        fullformula="formula",
                        x="matrix",
                        y="matrix",
                        terms="ANY",
                        xlevels="ANY",
                        contrasts="ANY",
                        preproc.x="function",
                        preproc.y="function",
                        "VIRTUAL"),
         prototype(formula=.~.,
                   fullformula=.~.,
                   preproc.x = function(x) x,
                   preproc.y = function(x) x))

## model-based clustering
setClass("FLXMC",
         representation(dist="character"),
         contains = "FLXM")

## regression
setClass("FLXMR",
         representation(offset="ANY"),
         contains = "FLXM")

setMethod("show", "FLXM",
function(object){
    cat("FlexMix model of type", object@name,"\n\nformula: ")
    print(object@formula)
    cat("Weighted likelihood possible:", object@weighted,"\n\n")
    if(nrow(object@x)>0){
        cat("Regressors:\n")
        print(summary(object@x))
    }
    if(nrow(object@y)>0){
        cat("Response:\n")
        print(summary(object@y))
    }
    cat("\n")
})

setClass("FLXcomponent",
         representation(df="numeric",
                        logLik="function",
                        parameters="list",
                        predict="function"))

setMethod("show", "FLXcomponent",
function(object){
    if(length(object@parameters)>0)
        print(object@parameters)
})
    


###**********************************************************

setClass("FLXP",
         representation(name="character",
                        formula="formula",
                        x="matrix",
                        fit="function",
                        refit="function",
                        coef="matrix",
                        df="function"),
         prototype(formula=~1, df = function(x, k, ...) (k-1)*ncol(x)))

setMethod("initialize", signature(.Object="FLXP"), function(.Object, ...) {
  .Object <- callNextMethod(.Object=.Object, ...)
  if (is.null(formals(.Object@refit))) .Object@refit <- .Object@fit
  .Object
})

setClass("FLXPmultinom",
         contains="FLXP")

setMethod("show", "FLXP",
function(object){
    cat("FlexMix concomitant model of type", object@name,"\n\nformula: ")
    print(object@formula)
    if(nrow(object@x)>0){
        cat("\nRegressors:\n")
        print(summary(object@x))
    }
    cat("\n")
})

###**********************************************************

setClass("FLXdist",
         representation(model="list",
                        prior="numeric",
                        components="list",
                        concomitant="FLXP",
                        formula="formula",
                        call="call",
                        k="integer"),
         validity=function(object) {
           (object@k == length(object@prior))
         },
         prototype(formula=.~.))

setClass("flexmix",
         representation(posterior="ANY",
                        weights="ANY",
                        iter="numeric",
                        cluster="integer",
                        logLik="numeric",
                        df="numeric",
                        control="FLXcontrol",
                        group="factor",
                        size="integer",
                        converged="logical",
                        k0="integer"),
         prototype(group=(factor(integer(0))),
                   formula=.~.),
         contains="FLXdist")

setMethod("show", "flexmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\nCluster sizes:\n")
    print(object@size)
    cat("\n")
    if(!object@converged) cat("no ")
    cat("convergence after", object@iter, "iterations\n")
})


###**********************************************************

setClass("summary.flexmix",
         representation(call="call",
                        AIC="numeric",
                        BIC="numeric",
                        logLik="logLik",
                        comptab="ANY"))

setMethod("show", "summary.flexmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    print(object@comptab, digits=3)
    cat("\n")
    print(object@logLik)
    cat("AIC:", object@AIC, "  BIC:", object@BIC, "\n")
    cat("\n")    
})

###**********************************************************

setClass("FLXMRglm",
         representation(family="character",
                        refit="function"),
         contains="FLXMR")

setClass("FLXR",
         representation(k="integer",
                        components = "list",
                        concomitant = "ANY",
                        call="call",
                        "VIRTUAL"))

setClass("FLXRoptim",
         representation(coef="vector",
                        vcov="matrix"),
         contains="FLXR")

setClass("FLXRmstep",
         contains="FLXR")
         
setMethod("show", signature(object = "FLXR"),
function(object) {
  cat("\nCall:", deparse(object@call,0.75*getOption("width")),
      sep="\n")
  cat("\nNumber of components:", object@k, "\n\n")
})

setMethod("summary", signature(object = "FLXRoptim"),
function(object, model = 1, which = c("model", "concomitant"), ...) {
  which <- match.arg(which)
  z <- if (which == "model") object@components[[model]] else object@concomitant
  show(z)
  invisible(object)
})

setMethod("summary", signature(object = "FLXRmstep"),
function(object, model = 1, which = c("model", "concomitant"), ...) {
  which <- match.arg(which)
  if (which == "model") {
    z <- object@components[[model]] 
    if (!is.null(z)) lapply(seq_along(z), function(k) {
      cat(paste("$", names(z)[k], "\n", sep = ""))
      printCoefmat(coef(summary(z[[k]])))
      cat("\n")
    })
  } else {
    z <- object@concomitant
    fitted.summary <- summary(z)
    k <- nrow(coef(fitted.summary)) + 1
    coefs <- lapply(2:k, function(n) {
      coef.p <- fitted.summary$coefficients[n - 1, , drop = FALSE]
      s.err <- fitted.summary$standard.errors[n - 1, , 
                                              drop = FALSE]
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
      coef.table <- t(rbind(coef.p, s.err, tvalue, pvalue))
      dimnames(coef.table) <- list(colnames(coef.p), c("Estimate", 
                                                       "Std. Error", "z value", "Pr(>|z|)"))
      new("Coefmat", coef.table)
    })
    names(coefs) <- paste("Comp", 2:k, sep = ".")
    print(coefs)
  }
  invisible(object)
})

setClass("Coefmat",
         contains = "matrix")

setMethod("show", signature(object="Coefmat"), function(object) {
  printCoefmat(object, signif.stars = getOption("show.signif.stars"))
})

###**********************************************************

setClass("FLXnested",
         representation(formula = "list",
                        k = "numeric"),
         validity = function(object) {
           length(object@formula) == length(object@k)
         })

setAs("numeric", "FLXnested",
      function(from, to) {
        new("FLXnested", formula = rep(list(~0), length(from)), k = from)
      })

setAs("list", "FLXnested",
      function(from, to) {
          z <- list2object(from, to)
      })

setAs("NULL", "FLXnested",
      function(from, to) {
        new(to)
      })

setMethod("initialize", "FLXnested", function(.Object, formula = list(), k = numeric(0), ...) {
    if (is(formula, "formula")) formula <- rep(list(formula), length(k))
    .Object <- callNextMethod(.Object, formula = formula, k = k, ...)
    .Object
})

###**********************************************************

setClass("FLXMRfix",
         representation(design = "matrix",
                        nestedformula = "FLXnested",
                        fixed = "formula",
                        segment = "matrix",
                        variance = "vector"),
         contains="FLXMR")

setClass("FLXMRglmfix",
         contains=c("FLXMRfix", "FLXMRglm"))

###**********************************************************

setClassUnion("listOrdata.frame", c("list", "data.frame"))

###**********************************************************

