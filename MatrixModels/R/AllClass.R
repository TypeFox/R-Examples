.onLoad <- function(lib, pkg) {
    if(is.null(getOption("max.print")))
	options(max.print = 10000)#-> show() of large matrices
}

## Model Matrix
setClass("modelMatrix",
	 representation(assign = "integer",
			contrasts = "list", "VIRTUAL"),
	 contains = "Matrix",
	 validity = function(object) {
	     if(length(object@assign) != (p <- ncol(object)))
		 return(gettextf("'%s' slot must be integer of length %d",
				 "assign", p))
	     contr <- object@contrasts
	     c.cl <- sapply(contr, class, USE.NAMES=FALSE)
	     if(length(nc <- names(contr)) != length(c.cl) || !all(nchar(nc) > 0))
		 return(gettextf("'%s' slot must be a correctly named list"))
	     ## TODO?  length(contrasts) < maximal value in 'assign' <= p
             contrCls <- c("character", "function", "matrix", "Matrix")
	     if(any(unlist(lapply(c.cl, function(cl) all(is.na(match(extends(cl), contrCls)))))))
		 return(gettextf("'%s' slot must be named list of contrast functions or their names, or matrices",
				 "contrasts"))
	     TRUE
	 })

setClass("sparseModelMatrix", representation("VIRTUAL"),
	 contains = c("CsparseMatrix", "modelMatrix"))
setClass("denseModelMatrix",  representation("VIRTUAL"),
	 contains = c("denseMatrix", "modelMatrix"))

## The currently only *actual* denseModelMatrix class:
setClass( "ddenseModelMatrix", contains = c("dgeMatrix", "ddenseMatrix", "denseModelMatrix"))
## here, add "ddense*": does not influence slots, but yields consistent superclass ordering

## The currently only *actual* sparseModelMatrix class:
setClass("dsparseModelMatrix", contains = c("dgCMatrix", "sparseModelMatrix"))

###------ Modules related to modelling --- basically in two parts --------------
###------ 1) "prediction-Module" -- currently in a sparse and dense flavor
###------ 2) "response-Module"

## Linear predictor modules, which consist of the model matrix, the
## coefficient vector and a triangular factor of the weighted model matrix.

## the super class contains the slots already;
setClass("predModule",
         representation(X = "modelMatrix", coef = "numeric", Vtr = "numeric",
                        fac = "CholeskyFactorization",
                        "VIRTUAL"))
## the sub classes specify more specific classes for the two non-trivial slots:
setClass("dPredModule", contains = "predModule",
	 representation(X = "ddenseModelMatrix", fac = "Cholesky"))
setClass("sPredModule", contains = "predModule",
	 representation(X = "dsparseModelMatrix", fac = "CHMfactor"))

## Response modules for models with a linear predictor, which can
## include linear models, generalized linear models, nonlinear models
## and generalized nonlinear models.

## y, offset and mu are as expected.  Note that length(offset) can be a multiple of length(y)
## weights are the prior weights
## sqrtrwt and sqrtXwt are the square roots of residual and X weights

setClass("respModule",
         representation(mu = "numeric",      # of length n
                        offset = "numeric",  # of length n * s
                        sqrtXwt = "matrix",  # of dim(.)  == (n, s)
                        sqrtrwt = "numeric", # sqrt(residual weights)
                        weights = "numeric", # prior weights
                        wtres = "numeric",
                        y = "numeric"),
         validity = function(object) {
             n <- length(object@y)
             if (any(n != sapply(lapply(c("weights","sqrtrwt","mu","wtres"
                     ), slot, object = object), length)))
                 return("lengths of weights, sqrtwt and mu must match length(y)")
             lo <- length(object@offset)
             if (!lo || lo %% n)
                 return("length(offset) must be a positive multiple of length(y)")
             if (length(object@sqrtXwt) != lo)
                 return("length(sqrtXwt) must equal length(offset)")
             if (nrow(object@sqrtXwt) != n)
                 return("nrow(sqrtXwt) != length(y)")
             TRUE
         })

setOldClass("family")

##' glm response module
setClass("glmRespMod",
         representation(family =  "family",
                        eta =    "numeric",
                        n =      "numeric"), # for evaluation of the aic
         contains = "respModule",
         validity = function(object) {
             if (length(object@eta) != length(object@y))
                 return("lengths of eta and y must match")
         })

##' nls response module
setClass("nlsRespMod",
         representation(nlenv = "environment",
                        nlmod = "call",
                        pnames = "character"),
         contains = "respModule",
         validity = function(object) {
             n <- length(object@y)
             N <- length(object@offset)
             s <- N %/% n
             lpn <- length(object@pnames)
             if (lpn != s) return(sprintf("length(pnames) = %d != s = %d", lpn, s))
             dd <- dim(object@sqrtXwt)
             if (!all(dd == c(n, s))) {
                 return(sprintf("dim(gradient) = (%d, %d), n = %d, s = %d",
                                dd[1], dd[2], n, s))
             }
             TRUE
         })

##' nglm response module
setClass("nglmRespMod", contains = c("glmRespMod", "nlsRespMod"))


### FIXME: move this eventually to 'methods':
##  -----
##' The mother class of all (S4 based)  (statistical / physical / ...) models in R:
setClass("Model", representation(call = "call", fitProps = "list",
                                 "VIRTUAL"))

##' Statistical models based on linear predictors
##' "glpModel" := General Linear Prediction Models
setClass("glpModel", representation(resp = "respModule", pred = "predModule"),
         contains = "Model")

rMod <-
    setRefClass("RespModule",
                fields = list(
                mu =      "numeric",    # of length n
                n =       "integer",    # for evaluation of the aic
                offset =  "numeric",    # of length n * s
                sqrtXwt = "matrix",     # of dim(.)  == (n, s)
                sqrtrwt = "numeric",    # sqrt(residual weights)
                weights = "numeric",    # prior weights
                wtres =   "numeric",
                y =       "numeric"),
                methods = list(
                initialize = function(...) {
                    initFields(...)
                    if (length(n) == 0L) n <<- length(y)
                    s <- 0L
		    ##currently fails at pkg INSTALL time: stopifnot(n > 0L)
                    if (length(weights) == 0L) weights <<- numeric(n) + 1
                    sqrtrwt <<- sqrt(weights)
                    if (any((dd <- dim(sqrtXwt)) < 1L))
                        sqrtXwt <<- matrix(1, ncol = 1L, nrow = n)
                    else {
                        stopifnot(nrow(sqrtXwt) == n)
                        s <- ncol(sqrtXwt)
                    }
                    swrk <- max(s, 1L)
                    if (length(offset) == 0) offset <<- numeric(n * swrk)
                    else {
                        so <- length(offset) %/% n
                        stopifnot(length(offset) %% n == 0, s == 0 || so == s)
                    }
                    wtres <<- mu <<- numeric(n) * NA_real_
                    .self
                },
                updateMu = function(gamma) {
                    gamma <- as.numeric(gamma)
                    stopifnot(length(gamma) == length(offset))
                    mu <<- gamma + offset
                    wtres <<- sqrtrwt * (y - mu)
                },
                updateWts = function() {}
                ))

rMod$lock("y","n","weights","offset")
glrMod <-
    setRefClass("GLMrespMod",
                fields = list(
                family =  "family",
                eta =    "numeric"),
                contains = "RespModule",
                methods = list(
                initialize = function(...) {
                    callSuper(...)
                    args <- list(...)
                    stopifnot("family" %in% names(args), is(args$family, "family"))
                    family <<- args$family
                    .self
                },
                updateMu = function(gamma) {
                    gamma <- as.numeric(gamma)
                    stopifnot(length(gamma) == length(offset))
                    mu <<- family$linkinv(eta <<- offset + gamma)
                    wtres <<- sqrtrwt * (y - mu)
                },
                updateWts = function() {
                    sqrtrwt <<- rtrwt <- sqrt(weights/family$variance(mu))
                    sqrtXwt[] <<- rtrwt * family$mu.eta(eta)
                    wtres <<- rtrwt * (y - mu)
                }
                ))
glrMod$lock("family")
