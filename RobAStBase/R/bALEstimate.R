###############################################################################
## Functions and methods for "ALEstimate" classes and subclasses
###############################################################################

setMethod("pIC", "ALEstimate", function(object) object@pIC)
setMethod("asbias", "ALEstimate", function(object) object@asbias)
setMethod("steps", "kStepEstimate", function(object) object@steps)
setMethod("Mroot", "MEstimate", function(object) object@Mroot)

setMethod("confint", signature(object="ALEstimate", method="missing"),
          function(object, method, level = 0.95) {
    objN <- paste(deparse(substitute(object)),sep="",collapse="")

    if(is.null(object@asvar)){ 
        cat(gettextf("Slot 'asvar' of object %s has not (yet) been filled.\n", objN))
        return(NULL) 
    }

    sd0 <- sqrt(diag(as.matrix(object@asvar))/object@samplesize)
    names(sd0) <- names(object@estimate)

### code borrowed from confint.default from package stats
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- .format.perc(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(object@estimate), 2),
                dimnames = list(names(object@estimate), pct))
    ci[] <- main(object) + sd0 %o% fac
### end of borrowed code

    new("Confint", type = gettext("asymptotic (LAN-based)"),
                   samplesize.estimate = object@samplesize,
                   call.estimate = object@estimate.call,
                   name.estimate = object@name,
                   trafo.estimate = object@trafo,
                   nuisance.estimate = nuisance(object),
                   fixed.estimate = fixed(object),
                   confint = ci)
})

setMethod("confint", signature(object="ALEstimate", method="symmetricBias"),
          function(object, method, level = 0.95) {
    objN <- paste(deparse(substitute(object)),sep="",collapse="")

    if(is.null(object@asvar)){ 
        cat(gettextf("Slot 'asvar' of object %s has not (yet) been filled.\n", objN))
        return(NULL) 
    }
    if(is.null(object@asbias)){ 
        cat(gettextf("Slot 'asbias' of object %s has not (yet) been filled.\n", objN))
        return(confint(object)) 
    }

    sd0 <- sqrt(diag(as.matrix(object@asvar))/object@samplesize)
    names(sd0) <- names(object@estimate)

### code borrowed from confint.default from package stats
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- .format.perc(a, 3)
    fac <- qnorm(a, mean = c(-object@asbias, object@asbias))
    ci <- array(NA, dim = c(length(object@estimate), 2),
                dimnames = list(names(object@estimate), pct))
    ci[] <- main(object) + sd0 %o% fac
### end of borrowed code

    new("Confint", type = c(
           gettext("asymptotic (LAN-based), uniform (bias-aware)\n"), 
           gettextf("for %s", name(method))
                           ),
                   samplesize.estimate = object@samplesize,
                   call.estimate = object@estimate.call,
                   name.estimate = object@name,
                   trafo.estimate = object@trafo,
                   nuisance.estimate = nuisance(object),
                   fixed.estimate = fixed(object),
                   confint = ci)
})

setMethod("confint", signature(object="ALEstimate", method="onesidedBias"),
          function(object, method, level = 0.95) {
    objN <- paste(deparse(substitute(object)),sep="",collapse="")

    if(is.null(object@asvar)){ 
        cat(gettextf("Slot 'asvar' of object %s has not (yet) been filled.\n", objN))
        return(NULL) 
    }
    if(is.null(object@asbias)){ 
        cat(gettextf("Slot 'asbias' of object %s has not (yet) been filled.\n", objN))
        return(confint(object)) 
    }

    sd0 <- sqrt(diag(as.matrix(object@asvar))/object@samplesize)
    names(sd0) <- names(object@estimate)

### code borrowed from confint.default from package stats
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- .format.perc(a, 3)
    if(method@sign == -1)
        M <- c(-object@asbias, 0)
    else
        M <- c(0, object@asbias)
    fac <- qnorm(a, mean = M)
    ci <- array(NA, dim = c(length(object@estimate), 2),
                dimnames = list(names(object@estimate), pct))
    ci[] <- main(object) + sd0 %o% fac
### end of borrowed code

    new("Confint", type = c(
           gettext("asymptotic (LAN-based), uniform (bias-aware)\n"), 
           gettextf("for %s", name(method))
                           ),
                   samplesize.estimate = object@samplesize,
                   call.estimate = object@estimate.call,
                   name.estimate = object@name,
                   trafo.estimate = object@trafo,
                   nuisance.estimate = nuisance(object),
                   fixed.estimate = fixed(object),
                   confint = ci)
})

setMethod("confint", signature(object="ALEstimate", method="asymmetricBias"),
          function(object, method, level = 0.95) {
    objN <- paste(deparse(substitute(object)),sep="",collapse="")

    if(is.null(object@asvar)){ 
        cat(gettextf("Slot 'asvar' of object %s has not (yet) been filled.\n", objN))
        return(NULL) 
    }
    if(is.null(object@asbias)){ 
        cat(gettextf("Slot 'asbias' of object %s has not (yet) been filled.\n", objN))
        return(confint(object)) 
    }

    sd0 <- sqrt(diag(as.matrix(object@asvar))/object@samplesize)
    names(sd0) <- names(object@estimate)

### code borrowed from confint.default from package stats
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- .format.perc(a, 3)
    fac <- qnorm(a, mean = c(-object@asbias, object@asbias)/method@nu)
    ci <- array(NA, dim = c(length(object@estimate), 2),
                dimnames = list(names(object@estimate), pct))
    ci[] <- main(object) + sd0 %o% fac
### end of borrowed code

    nuround <- round(nu,3)
    new("Confint", type = c(
           gettext("asymptotic (LAN-based), uniform (bias-aware)\n"), 
           gettextf("for %s with nu =(%f,%f)", 
                     name(method), nuround[1], nuround[2])
                           ),
                   samplesize.estimate = object@samplesize,
                   call.estimate = object@estimate.call,
                   name.estimate = object@name,
                   trafo.estimate = object@trafo,
                   nuisance.estimate = nuisance(object),
                   fixed.estimate = fixed(object),
                   confint = ci)
})
