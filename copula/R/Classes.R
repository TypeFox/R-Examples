## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### basic copula class #########################################################


setClass("copula",
         representation(dimension = "integer", # as for "nacopula"
                        parameters = "numeric",
                        param.names = "character",
                        param.lowbnd = "numeric",
                        param.upbnd = "numeric",
                        ## TODO: "a vector" of intervals" paraInterval = "maybeInterval", # [.,.]  (.,.], etc .. parameter interval
                        fullname = "character",
                        "VIRTUAL"),
         prototype = prototype(dimension = 2L, parameters = NA_real_),
         validity = ##' Check validity of "copula"
         function(object) {
	     dim <- object@dimension # "integer" by definition
	     if (length(dim) != 1L) return("'dim' must be an integer (>= 2)")
	     if (dim < 2) return("dim must be >= 2")
             param <- object@parameters
             upper <- object@param.upbnd
             lower <- object@param.lowbnd
             lp <- length(param)
             if (lp != length(upper) && length(upper) != 1)
                 return("Parameter and upper bound have non-equal length")
             if (lp != length(lower) && length(lower) != 1)
                 return("Parameter and lower bound have non-equal length")
             intervChk <- ## TODO: mkParaConstr(object@paraInterval)
                 function(par) all(is.na(param) | (lower <= param & param <= upper))
             ina.p <- is.na(param)
             if(!all(ina.p)) {
		 ##if(any(ina.p)) return("some (but not all) parameter values are  NA")
                 if(!intervChk(param)) return("Parameter value(s) out of bound")
             }

	     ## want to allow (all) NA parameters:
	     TRUE
         })

## general methods for copula
setGeneric("dCopula", function(u, copula, log=FALSE, ...) {
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    u.is.out <- outside.01(u, strictly=FALSE)## on.boundary _or_ outside
    if(any.out <- any(u.is.out, na.rm=TRUE))
	u[] <- pmax(0, pmin(1, u)) # <- "needed", as some methods give error
    r <- standardGeneric("dCopula")
    if(any.out) ## on boundary _or_ outside cube  ==> zero mass :
	r[u.is.out & !is.na(u.is.out)] <- if(log) -Inf else 0.
    r
})
setGeneric("pCopula", function(u, copula, ...) {
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    ## here as well, 'outside' and 'on-boundary' are equivalent:
    u[] <- pmax(0, pmin(1, u))
    standardGeneric("pCopula")
})
setGeneric("rCopula", function(n, copula, ...) standardGeneric("rCopula"))
setGeneric("tau", function(copula, ...) standardGeneric("tau"))
setGeneric("rho", function(copula, ...) standardGeneric("rho"))
setGeneric("tailIndex", function(copula, ...) standardGeneric("tailIndex"))
setGeneric("iTau", function(copula, tau, ...) standardGeneric("iTau"))
setGeneric("iRho", function(copula, rho, ...) standardGeneric("iRho"))
setGeneric("dTau", function(copula, ...) standardGeneric("dTau"))
setGeneric("dRho", function(copula, ...) standardGeneric("dRho"))

setGeneric("dTauFun", function(copula) standardGeneric("dTauFun"))
setGeneric("dRhoFun", function(copula) standardGeneric("dRhoFun"))

## Deprecated:
calibKendallsTau <- function(copula, tau) { .Deprecated("iTau"); iTau(copula,tau) }
calibSpearmansRho <- function(copula, rho) { .Deprecated("iRho"); iRho(copula,rho) }
kendallsTau <- function(copula) { .Deprecated("tau"); tau(copula) }
spearmansRho <- function(copula) { .Deprecated("rho"); rho(copula) }
genInv <- function(copula, s) { .Deprecated("psi"); psi(copula,s) }
genFun <- function(copula, u) { .Deprecated("iPsi"); iPsi(copula, u) }
genFunDer1 <- function(copula, u){ .Deprecated("diPsi"); diPsi(copula, u) }
genFunDer2 <- function(copula, u){ .Deprecated("diPsi(*, degree=2)"); diPsi(copula, u, degree=2) }

AfunDer <- function(copula, w) { .Deprecated("dAdu"); dAdu(copula, w) }
Afun    <- function(copula, w) { .Deprecated("A"); A(copula, w) }

pcopula <- function(copula, u, ...) { .Deprecated("pCopula"); pCopula(u, copula) }
dcopula <- function(copula, u, ...) { .Deprecated("dCopula"); dCopula(u, copula, ...) }
rcopula <- function(copula, n, ...) { .Deprecated("rCopula"); rCopula(n, copula, ...) }


### elliptical copulas, contains normalCopula and tCopula ######################

## NOTE: This is related to  npar.ellip() in ./ellipCopula.R
validRho <- function(dispstr, dim, lenRho) {
    switch(dispstr, ## checking for correct 'dispstr'
	   "ar1" =, "ex" = {
	       if (lenRho != 1)
		   return(gettextf("'rho' parameter should have length 1 for 'dispstr' = \"%s\"",
				   dispstr))
	   },
	   "un" = {
	       if (lenRho != (L <- dim * (dim - 1) / 2))
		   return(gettextf("'rho' parameter should have length dim * (dim - 1) / 2 (= %d) for 'dispstr' = \"%s\"",
				   L, dispstr))
	   },
	   "toep" = {
	       if (lenRho != dim - 1)
		   return(gettextf("'rho' parameter should have length dim-1 (= %d) for 'dispstr' = \"%s\"",
				   dim-1, dispstr))
	   },
	   ## otherwise
	   return("'dispstr' not supported (yet)"))

    TRUE
}

validEllipCopula <- function(object) {
  rho <- object@getRho(object)
  validRho(dispstr=object@dispstr, dim=object@dimension,
           length(rho))
}

setClass("ellipCopula", contains = "copula",
	 representation(dispstr = "character", getRho="function", "VIRTUAL"),
         validity = validEllipCopula)

if(FALSE) # not yet needed -- validEllipCopula() is used anyway
##' normal copula
validNormalCopula <- function(object) {
  ## can do more if needed here
}
setClass("normalCopula", contains = "ellipCopula")
         ## validity = validNormalCopula)

if(FALSE)## not needed
## t copula
validTCopula <- function(object) {
  ## df inside boundaries is checked in "copula" validity
}

setClass("tCopula", representation(df = "numeric", df.fixed = "logical"),
         contains = "ellipCopula"
         ## , validity = validTCopula
         )


## methods for ellipCopula??

### Archimedean copulas, contains AMH, Clayton, Frank, Gumbel, ... #############

setClass("archmCopula", representation(exprdist = "expression", "VIRTUAL"),
	 contains = "copula")

## clayton copula
setClass("claytonCopula", contains = "archmCopula")

## gumbel copula, also an ev copula

## frank copula
setClass("frankCopula", contains = "archmCopula")

## amh copula
setClass("amhCopula", contains = "archmCopula")

## Joe copula
setClass("joeCopula", contains = "archmCopula")

## methods for archmCopulas
setGeneric("psi", function(copula, s) standardGeneric("psi"))
##FIXME 'log' compulsory:
##setGeneric("iPsi", function(copula, u, log, ...) standardGeneric("iPsi"))
setGeneric("iPsi", function(copula, u, ...) standardGeneric("iPsi"))
setGeneric("dPsi", function(copula, s, ...) standardGeneric("dPsi"))
setGeneric("diPsi", function(copula, u, degree=1, log=FALSE, ...) standardGeneric("diPsi"))


### Extreme value copulas, contains galambos, husler-reiss, gumbel, ... ########

setClass("evCopula", representation("VIRTUAL"), contains = "copula")

## galambos copula
setClass("galambosCopula", representation(exprdist = "expression"),
         contains = "evCopula")

## gumbel copula, also an archm copula;
setClass("gumbelCopula", contains = list("archmCopula", "evCopula"))

## husler-reiss copula
setClass("huslerReissCopula",representation(exprdist = "expression"),
         contains = "evCopula")

## tawn copula; does not offer full range of dependence
setClass("tawnCopula", representation(exprdist = "expression"),
         contains = "evCopula")

## tEV copula
setClass("tevCopula", representation(df = "numeric", df.fixed = "logical"),
         contains = "evCopula")

setGeneric("A", function(copula, w) standardGeneric("A"))
setGeneric("dAdu", function(copula, w) standardGeneric("dAdu"))
setGeneric("dAdtheta", function(copula, w) standardGeneric("dAdtheta"))

### independent copula class ###################################################

## it should contain all three *virtual* superclasses,
## but we don't want it to inherit "funny slots"
setClass("indepCopula", contains = c("evCopula", "archmCopula"))


### Other copulas ##############################################################

## Farlie-Gumbel-Morgenstern multivariate copula
setClass("fgmCopula", representation(exprdist = "expression"),
         contains = "copula",
         ## verify that the pdf is positive at each vertex of [0,1]^dim
         validity = function(object) {
             dim <- object@dimension
             if (dim == 2)
                 return(TRUE)
             param <- object@parameters
             valid <- .C(validity_fgm,
                         as.integer(dim),
                         as.double(c(rep(0,dim+1),param)),
                         valid = integer(1))$valid
             if (valid == 0)
                 return("Bad vector of parameters")
             else
                 return(TRUE)
         })


## plackett copula
setClass("plackettCopula",representation(exprdist = "expression"),
         contains = "copula")

### Multivariate distibution via copula ########################################

validMvdc <- function(object) {
    dim <- object@copula@dimension
    if(!is.finite(dim) || dim < 2)
	return("'dimension' must be integer >= 2")
    if(dim != length(object@margins))
	return("'dimension' does not match margins' length")
    if(dim != length(pm <- object@paramMargins))
	return("'dimension' does not match paraMargins' length")
    if(!all(vapply(pm, function(e) is.list(e) || is.numeric(e), NA)))
	return("'paramMargins' elements must all be list()s or numeric vectors")
    okNms <- function(nms) !is.null(nms) && all(nzchar(nms))
    if(object@marginsIdentical) {
	if(!all(object@margins[1] == object@margins[-1]))
	    return("margins are not identical")
	pm1 <- pm[[1]]
	for(i in 2:dim) {
	    if(!identical(pm1, pm[[i]]))
		return("margins are not identical")
	}
	if(length(pm1) > 0 && !okNms(names(pm1)))
	    return("'paramMargins' must be named properly")
    }
    else ## not identical margins: check each
	for(i in seq_len(dim)) {
	    pmi <- pm[[i]]
	    if(length(pmi) > 0 && !okNms(names(pmi)))
		return(gettextf("'paramMargins[[%d]]' must be named properly", i))
	    ## TODO(?): check more similar to (/ instead of) those in mvdc() --> ./mvdc.R
	}
    TRUE
}## validMvdc()

setClass("mvdc",
	 representation(copula = "copula",
			margins = "character",
			paramMargins = "list",
			marginsIdentical = "logical"),
	 validity = validMvdc)

## methods like {dpr}mvdc are defined in mvdc.R

## A fitted multivariate distribution -- "generic mother class",
## "fitCopula" and "fitMvdc" will inherit from it:
setClass("fittedMV",
	 representation(estimate = "numeric",
			var.est = "matrix", ## FIXME 'vcov'
			loglik = "numeric",
			nsample = "integer",
			method = "character",
			## convergence = "integer",
			fitting.stats = "list"))

setGeneric("paramNames", function(x) standardGeneric("paramNames"))
## paramNames() methods: provided separately for "fitCopula", "fitMvdc"

coef.fittedMV <- function(object, ...)
    setNames(object@estimate, paramNames(object))

nobs.fittedMV <- function(object, ...) object@nsample

vcov.fittedMV <- function(object, ...) {
    pNms <- paramNames(object)
    structure(object@var.est, dimnames = list(pNms, pNms))
}

logLik.fittedMV <- function(object, ...) {
    val <- object@loglik
    attr(val, "nobs") <- object@nsample
    attr(val, "df") <- length(object@estimate)
    class(val) <- "logLik"
    val
}

###-------------------------- Glue   "copula" <-> "nacopula"

##' The mother of all copula classes:
setClassUnion("Copula",
              members = c("copula", "nacopula"))
## NB: "acopula" *not* : It has no dimension, is rather a family object

##' does the copula have 'df' as parameter?
has.par.df <- function(cop, classDef = getClass(class(cop)),
                       isEllip = extends(classDef, "ellipCopula")) {
  ((isEllip && extends(classDef, "tCopula")) ||
   extends(classDef, "tevCopula")) && !cop@df.fixed
}
