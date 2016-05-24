###############################################################################
## confint methods
###############################################################################

setMethod("confint", signature(object="ANY", method="missing"),
        function(object, method, parm, level = 0.95, ...) {
        if(hasArg(parm))
           stats::confint(object = object, parm = parm, level = level, ...)
        else
           stats::confint(object = object, level = level, ...)
})

setMethod("confint", signature(object="Estimate", method="missing"),
          function(object, method, level = 0.95) {
   objN <- paste(deparse(substitute(object)),sep="",collapse="")

   asm <- asvar(object)
   if(is.null(asm)){
      cat(gettextf("Slot 'asvar' of object %s has not (yet) been filled.\n",
          objN))
      return(NULL)
   }
   sd0 <- sqrt(diag(asm)/object@samplesize)
   names(sd0) <- names(object@estimate)

### code borrowed from confint.default from package stats
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- stats:::format.perc(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(object@estimate), 2),
                dimnames = list(names(object@estimate), pct)
                )
    ci[] <- main(object) + sd0 %o% fac
### end of borrowed code
    new("Confint", type = gettext("asymptotic (CLT-based)"),
                   samplesize.estimate = object@samplesize,
                   call.estimate = object@estimate.call,
                   name.estimate = object@name,
                   trafo.estimate = object@trafo,
                   nuisance.estimate = nuisance(object),
                   fixed.estimate = fixed(object),
                   confint = ci)
})

if(require(stats4)){
setMethod("confint", signature(object="mle", method="missing"),
        function(object, method, parm, level = 0.95, ...) {
        if(hasArg(parm))
           stats4::confint(object = object, parm = parm, level = level, ...)
        else
           stats4::confint(object = object, level = level, ...)
})
setMethod("confint", signature(object="profile.mle", method="missing"),
        function(object, method, parm, level = 0.95, ...) {
        if(hasArg(parm))
           stats4::confint(object = object, parm = parm, level = level, ...)
        else
           stats4::confint(object = object, level = level, ...)
})
}
