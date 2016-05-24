getLightFit <- function(#Makes Lighter Fitted Object
### Makes a lighter version of a fitted object by removing elements containing data.     
    fit
### A fitted object.
    ) {
  ##details<< Most of  the space used by a fitted object  is not necessary for
  ##prediction.   This  concerns, for  instance,  the "residuals",  "effects",
  ##"fitted.values",  and  "model"  entries   of  a  linear  model  fitted  by
  ##\code{lm}.   These entries  can thus  be removed  from the  object without
  ##affecting  the  model  predictions.    This  function  is  currently  only
  ##implemented for fitted objects that  derive from class 'lm' or 'rpart'. It
  ##is mainly for internal use.

  knownClasses <- c("lm", "glm", "rpart")
  known <- sapply(knownClasses, FUN=function(cls) inherits(fit, cls))
  if (!any(known)) {
    throw("Object type not supported: ", class(fit))
  }
  elts <- NULL
  if (known[["lm"]]) {
    elts <- c(elts, "residuals", "effects", "fitted.values", "model")
    fit$qr$qr <- NA
  }
  if (known[["glm"]]) {     
    elts <- c(elts, "linear.predictors", "weights", "prior.weights", "y", "data")
  }
  if (known[["rpart"]]) {     
    elts <- c(elts, "where", "y")
    fit
  }
  fit[elts] <- NA
  fit
### Returns the same object as the input without the entries that contain data.
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

