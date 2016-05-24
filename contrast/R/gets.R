getCoefficients <- function(fit, ...) UseMethod("getCoefficients")

# the default version of the getCoefficients generic function simply
# calls coef, and sanity checks the results.
getCoefficients.default <- function(fit, ...)
{
   tmp <- coef(fit)

   # sanity check the return value.  if this fails, then we probably need a
   # special version of getCoefficients for this class of fit objects.
   if (! is.vector(tmp) || is.list(tmp) || is.matrix(tmp))
      stop(paste('coef function returned bad type for fit object of class:',
                 paste(class(fit), collapse=" ")))
   tmp
}

getCoefficients.lme <- function(fit, ...)
{
  library(nlme)
   fixef(fit)
}
