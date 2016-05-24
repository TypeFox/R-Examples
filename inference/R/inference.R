##' Extract inferential information from different statistical models.
##'
##' This package provides functions to extract point estimates, standard
##' errors, confidence intervals, p-values, and sample size from a fitted
##' model in a matrix-like object.  The purpose is to have all inferential
##' numbers be readily accessible, especially in the construction
##' of summary tables (R -> LaTeX -> html -> Word) for publication and
##' collaboration.
##'
##' @name inference-package
##' @docType package
##' @author Vinh Nguyen \email{vinhdizzo at gmail dot com}
##' @keywords statistical models point estimates confidence intervals
##' p-values sample size inference
##' @import methods
##' @examples
##' infer(lm(rnorm(100) ~ runif(100)))
NULL
## don't use @rdname alone in above.  will just use @name


##' An S4 class that stores inferential values of a fitted model object.
##'
##' An S4 class that inherits from the \code{\link[=structure-class]{matrix}}
##' class in the \code{methods} package; see \code{class?matrix}.  Rows
##' correspond to different coefficients and columns consist of point
##' estimates (point.est), confidence intervals (ci.lo and ci.hi),
##' p-values (p.value), and sample size (n).
##'
##' @rdname inference-class
##' @docType class
##' @slot .Data
##' @slot model String specifying class of model fit, such as "lm".
##' @slot sample.size Sample size used in model fit.
##' @slot robust.se Boolean indicator whether robust standard errors were
##' used.
##' @slot two.sided Boolean indicator whether p-values corresond to a
##' two-sided test or one-sided.
##' @slot ci.level Confidence level, e.g., 0.95.
##' @slot scale Scale of point estimates; defaults to "beta".
##' @slot others List containing other information about the model;
##' eg, summary of cluster size for \code{\link[gee]{gee}} and \code{\link[nlme]{lme}} objects;
##' number of events for \code{\link[survival]{coxph}} objects.
##' @exportClass inference
setClass(Class="inference"
         , representation=representation(model="character"
             , sample.size="numeric"
             , robust.se="logical"
             , two.sided="logical"
             , ci.level="numeric"
             , scale="character"
             , others="list")
         , contains=c("matrix"))

##' Show/print \code{\linkS4class{inference}} object.
##'
##' \code{\link[methods]{show}} method for objects made using the \code{\link{infer}} function.
##' @rdname show,inference-method
##' @aliases print.inference
##' @param object \code{\linkS4class{inference}} object.
##' @return Nothing.
##' @author Vinh Nguyen
setMethod("show", "inference", function(object){ print(slot(object, ".Data"))})

##' Inference for fitted model objects.
##'
##' Extract point estimates, standard errors, confidence intervals,
##' p-values, and sample size.
##' @rdname infer,-methods
##' @aliases infer infer,-method infer,lm-method infer,glm-method infer,coxph-method infer,gee-method infer,lme-method infer,mer-method
##' @docType methods
##' @usage infer(fitobj, vars, robust.se=TRUE, two.sided=TRUE, ci.level=0.95, ...)
##' @param fitobj Fitted model object, such as those of class \code{\link[stats]{lm}}.
##' @param vars Vector of variable names to obtain inference information
##' for.  If not specified, all variables in the fitted model will be used.
##' in the fitted model.
##' @param robust.se Boolean indicator for whether robust standard
##' errors should be use.  Defaults to \code{TRUE}.
##' @param two.sided Boolean indicator for whether p-values should
##' correspond to a two-sided test or one-sided.  Defaults to
##' \code{TRUE}.
##' @param ci.level Confidence level.  Defaults to 0.95.
##' @param ... Not used.
##' @return S4 \code{\linkS4class{inference}} object.
##' @examples
##' infer(lm(rnorm(100) ~ runif(100)))
##' @exportMethod infer
##' @author Vinh Nguyen
setGeneric(name="infer", function(fitobj, vars, robust.se=TRUE, two.sided=TRUE, ci.level=0.95, ...) standardGeneric("infer"))
##setGeneric(name="infer", function(fitobj, ...) standardGeneric("infer"))

## add infer,lm-method to @aliases in infer documentation
##' @nord
setMethod("infer", signature(fitobj="lm"), function(fitobj, vars, robust.se=TRUE, two.sided=TRUE, ci.level=0.95)
{
  if(missing(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    require(sandwich)
    se <- sqrt(diag(sandwich(fitobj))[vars])
  } else{
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(fitobj$residuals)
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, scale="beta", others=list("empty"))
  return(rslt)
})

## add infer,glm-method to @aliases in infer documentation
##' @nord
setMethod("infer", signature(fitobj="glm"), function(fitobj, vars, robust.se=TRUE, two.sided=TRUE, ci.level=0.95)
{
  if(missing(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    require(sandwich)
    se <- sqrt(diag(sandwich(fitobj))[vars])
  } else{
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(fitobj$residuals)
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, scale="beta", others=list("empty"))
  return(rslt)
})

## add infer,coxph-method to @aliases in infer documentation
##' @nord
setMethod("infer", signature(fitobj="coxph"), function(fitobj, vars, robust.se=TRUE, two.sided=TRUE, ci.level=0.95)
{
  if(missing(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    require(sandwich)
    if(any(names(fitobj) == "naive.var")){
      se <- sqrt(diag(vcov(fitobj))[vars])
    } else{
      se <- sqrt(diag(sandwich(fitobj))[vars])
    }
  } else{
    if(any(names(fitobj) == "naive.var")){
      warning("Only robust standard errors available due to robust=TRUE specification in coxph.")
    }
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(fitobj$residuals)
  n.events <- sum(fitobj$y[, 2])
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, scale="beta", others=list(n.events=n.events))
  return(rslt)
})


## add infer,gee-method to @aliases in infer documentation
##' @nord
setMethod("infer", signature(fitobj="gee"), function(fitobj, vars, robust.se=TRUE, two.sided=TRUE, ci.level=0.95)
{
  if(missing(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    ##require(sandwich)
    se <- sqrt(diag(fitobj$robust.variance)[vars])
  } else{
    se <- sqrt(diag(fitobj$naive.variance)[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(unique(fitobj$id))
  nObs <- fitobj$nobs
  summaryClusters <- summary(tapply(fitobj$id, fitobj$id, length))
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, scale="beta", others=list(nObs=nObs, summaryClusters=summaryClusters))
  return(rslt)
})

## add infer,lme-method to @aliases in infer documentation
##' @nord
setMethod("infer", signature(fitobj="lme"), function(fitobj, vars, robust.se=FALSE, two.sided=TRUE, ci.level=0.95)
{
  if(missing(vars)) vars <- names(coef(fitobj))
  point.est <- fixed.effects(fitobj)[vars] ##coef(fitobj)[vars]
  if(robust.se){
    ##require(sandwich)
    ##se <- sqrt(diag(sandwich(fitobj))[vars])
    stop("Robust standard errors are not available with Linear Mixed Effects Models.")
  } else{
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- fitobj$dims$ngrps["Subject"]
  nObs <- fitobj$dims$N
  summaryClusters <- summary(tapply(fitobj$groups[, 1], fitobj$groups[, 1], length))
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, scale="beta", others=list(nObs=nObs, summaryClusters=summaryClusters))
  return(rslt)
})

## add infer,mer-method to @aliases in infer documentation
##' @nord
setMethod("infer", signature(fitobj="mer"), function(fitobj, vars, robust.se=FALSE, two.sided=TRUE, ci.level=0.95)
{
  if(missing(vars)) vars <- colnames(coef(fitobj)[[1]])
  point.est <- fitobj@fixef[vars] ##coef(fitobj)[vars]
  SEfitobj <- sqrt(diag(vcov(fitobj)))
  names(SEfitobj) <- names(fitobj@fixef)
  if(robust.se){
    ##require(sandwich)
    ##se <- sqrt(diag(sandwich(fitobj))[vars])
    stop("Robust standard errors are not available with Linear Mixed Effects Models.")
  } else{
    se <- se <- SEfitobj[vars]
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- slot(fitobj, "dims")["q"]
  nObs <- slot(fitobj, "dims")["n"]
  summaryClusters <- summary(tapply(slot(fm1, "flist")[, 1], slot(fm1, "flist")[, 1], length))
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, scale="beta", others=list(nObs=nObs, summaryClusters=summaryClusters))
  return(rslt)
})



##' \code{\link[base]{transform}} method for class inference
##'
##' Transform the point estimates, confidence intervals, and standard
##' errors based on the delta method.  This builds on the S3 generic
##' function \code{\link[base]{transform}} from the \code{base} package.
##'
##' It can be used to get the hazard ratio scale in inference objects
##' created from \code{\link[survival]{coxph}} objects and the odds ratio scale from
##' logistic regression created from \code{\link[stats]{glm}}
##' (both using \code{f=exp, f.prime=exp}).
##' @title Transformation of point estimates
##' @rdname transform.inference
##' @usage transform.inference(`_data`, f, f.prime, ...)
##' @param _data Object of class \code{\linkS4class{inference}}.
##' @param f Function to transform the point estimates and confidence intervals; e.g., \code{\link[base]{exp}}.
##' @param f.prime Derivative of \code{f} in order to compute the standard
##' error of the transformed point estimates based on the delta method.
##' @param ... Nothing.
##' @return Object of class \code{\linkS4class{inference}}.
##' @author Vinh Nguyen
transform.inference <- function(`_data`, f, f.prime, ...)
{
  infObj <- infObjTrans <- `_data`
  stopifnot(is.function(f), class(infObj)=="inference")
  infObjTrans[, "point.est"] <- f(infObj[, "point.est"])
  infObjTrans[, "ci.lo"] <- f(infObj[, "ci.lo"])
  infObjTrans[, "ci.hi"] <- f(infObj[, "ci.hi"])
  if(missing(f.prime) | is.null(f.prime)){
    warning("f.prime not provided.  Standard errors for the resulting 'inference' object does NOT correspond to that of the transformed point estimates.")
  } else{
    infObjTrans[, "se"] <- f.prime(infObj[, "point.est"]) * infObj[, "se"]
  }
  scale <- paste(as.character(substitute(f)), "( ", slot(infObj, "scale"), " )", sep="")
  slot(infObjTrans, "scale") <- scale
  return(infObjTrans)
}

## infObj <- infer(lm(runif(100) ~ rnorm(100)))
## library(roxygen)
## package.skeleton("inference2", code_files="inference.R", force=TRUE)
## roxygenize(package.dir="inference", roxygen.dir="inference", copy.package=FALSE, unlink.target=FALSE)
## roxygenize(package.dir="inference", roxygen.dir="inference", copy.package=FALSE, unlink.target=FALSE, use.Rd2=TRUE) ## will try to document all function declarations, even if I did not document it using roxygen
## roxygenize(package.dir="pkg", roxygen.dir="pkg", copy.package=FALSE, unlink.target=FALSE, use.Rd2=TRUE) ## use this for S4
## R CMD roxygen -d -s pkg
