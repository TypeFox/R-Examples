##
## Generic definition for estimating GO-GARCH models
##
setGeneric("goest", function(object, ...) standardGeneric("goest"))
##
## Generic definition for extracting Euler angles
##
setGeneric("angles", function(object, ...) standardGeneric("angles"))
##
## Generic definition for extracting the conditional variances
##
setGeneric("cvar", function(object, ...) standardGeneric("cvar"))
##
## Generic definition for extracting the conditional covariances
##
setGeneric("ccov", function(object, ...) standardGeneric("ccov"))
##
## Generic definition for extracting the conditional correlations
##
setGeneric("ccor", function(object, ...) standardGeneric("ccor"))
##
## Generic definition for extracting convergence codes
##
setGeneric("converged", function(object, ...) standardGeneric("converged"))
##
## Generic definition for extracting object@M for objects of class Orthom
##
setGeneric("M", function(object, ...) standardGeneric("M"))
##
## Setting Generics for coef, formula, logLik,
## plot, predict, residuals, resid, summary, t
## and update
##
setGeneric("coef")
setGeneric("formula")
setGeneric("logLik")
setGeneric("plot")
setGeneric("predict")
setGeneric("residuals")
setGeneric("resid")
setGeneric("summary")
setGeneric("t")
setGeneric("update")
