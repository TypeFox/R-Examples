##*****************************************************************
## Some methods for covariance classes.
##
## XXX plot methods not complete
##
##*****************************************************************



##=================================================================
## summary should in the future return objects with the right class
##=================================================================

if(!isGeneric("summary")) {
  setGeneric(name = "summary",
             def = function(object, ...) standardGeneric("summary"))
}

summary.covIso <- function(object, ...) {
  ans <- list()
  ans[["input"]] <- object@d
  ans[["inputnames"]] <- object@var.names
  ans[["kernelname"]] <- object@name
  ans 
  
}

setMethod("summary", signature(object = "covIso"),
          definition = summary.covIso
          )

summary.covTensorProduct <- function(object, ...) {
  ans <- list()
  ans[["input"]] <- object@d
  ans[["inputnames"]] <- object@var.names
  ans[["kernelname"]] <- object@name
  ans[["ranges"]] <- object@range.val
  names(ans[["ranges"]]) <- object@var.names
  ans 
  
}

setMethod("summary", signature(object = "covTensorProduct"),
          definition = summary.covTensorProduct
          )

summary.covAffineScaling <- function(object, ...) {
  ans <- list()
  ans[["input"]] <- object@d
  ans[["inputnames"]] <- object@var.names
  ans[["kernelname"]] <- object@name
  ans[["knots"]] <- object@knots
  ans[["eta"]] <- object@eta
  names(ans[["eta"]]) <- object@var.names
  ans 
}

setMethod("summary", signature(object = "covAffineScaling"),
          definition = summary.covAffineScaling
          )

setMethod("summary", signature(object = "covScaling"),
          definition = summary.covAffineScaling
          )


##=================================================================
## ninput
##=================================================================

if(!isGeneric("ninput")) {
  setGeneric(name = "ninput",
             def = function(x) standardGeneric("ninput"))
}

setMethod("ninput", signature(x = "covIso"), function(x) x@d)
setMethod("ninput", signature(x = "covTensorProduct"), function(x) x@d)
setMethod("ninput", signature(x = "covAffineScaling"), function(x) x@d)
setMethod("ninput", signature(x = "covScaling"), function(x) x@d)


##=================================================================
## inputnames
##=================================================================

if(!isGeneric("inputnames")) {
  setGeneric(name = "inputnames",
             def = function(x) standardGeneric("inputnames"))
}

setMethod("inputnames", signature(x = "covIso"), function(x) x@var.names)
setMethod("inputnames", signature(x = "covTensorProduct"), function(x) x@var.names)
setMethod("inputnames", signature(x = "covAffineScaling"), function(x) x@var.names)
setMethod("inputnames", signature(x = "covScaling"), function(x) x@var.names)

# below: beta version. One must also replace the names in X and y, for consistency.
# 
# setGeneric("inputnames<-", function(x, value){ standardGeneric("inputnames<-") })
# 
# setReplaceMethod("inputnames", signature(x = "covTensorProduct", value = "character"),
#                  function(x, value) {
#                    x@var.names <- value
#                    return(x)}
#                  )
# setReplaceMethod("inputnames", signature(x = "covIso", value = "character"),
#                  function(x, value) {
#                    x@var.names <- value
#                    return(x)}
#                  )
# setReplaceMethod("inputnames", signature(x = "covAffineScaling", value = "character"),
#                  function(x, value) {
#                    x@var.names <- value
#                    return(x)}
#                  )

##=================================================================
## kernelname
##=================================================================

if(!isGeneric("kernelname")) {
  setGeneric(name = "kernelname",
             def = function(x) standardGeneric("kernelname"))
}

setMethod("kernelname", signature(x = "covIso"), function(x) x@name)
setMethod("kernelname", signature(x = "covTensorProduct"), function(x) x@name)
setMethod("kernelname", signature(x = "covAffineScaling"), function(x) x@name)
setMethod("kernelname", signature(x = "covScaling"), function(x) x@name)

##=================================================================
## nuggetvalue
##=================================================================

if(!isGeneric("nuggetvalue")) {
  setGeneric(name = "nuggetvalue",
             def = function(x) standardGeneric("nuggetvalue"))
}

setMethod("nuggetvalue", signature(x = "covIso"), function(x) x@nugget)
setMethod("nuggetvalue", signature(x = "covTensorProduct"), function(x) x@nugget)
setMethod("nuggetvalue", signature(x = "covAffineScaling"), function(x) x@nugget)
setMethod("nuggetvalue", signature(x = "covScaling"), function(x) x@nugget)
setMethod("nuggetvalue", signature(x = "covUser"), function(x) x@nugget)


setGeneric("nuggetvalue<-",function(x, value){ standardGeneric("nuggetvalue<-") })

nuggetvalueFun <- function(x, value) {
  if (length(value)>0) {
    x@nugget <- value
    x@nugget.flag <- TRUE
  } else {
    x@nugget <- numeric(0)
    x@nugget.flag <- FALSE
  }
  validObject(x)
  return(x)
}

setReplaceMethod("nuggetvalue", signature(x = "covTensorProduct", value = "numeric"), nuggetvalueFun)
setReplaceMethod("nuggetvalue", signature(x = "covIso", value = "numeric"), nuggetvalueFun)             
setReplaceMethod("nuggetvalue", signature(x = "covAffineScaling", value = "numeric"), nuggetvalueFun)
setReplaceMethod("nuggetvalue", signature(x = "covScaling", value = "numeric"), nuggetvalueFun)
setReplaceMethod("nuggetvalue", signature(x = "covUser", value = "numeric"), nuggetvalueFun)

##=================================================================
## nuggetflag
##=================================================================

if(!isGeneric("nuggetflag")) {
  setGeneric(name = "nuggetflag",
             def = function(x) standardGeneric("nuggetflag"))
}

setMethod("nuggetflag", signature(x = "covIso"), function(x) x@nugget.flag)
setMethod("nuggetflag", signature(x = "covTensorProduct"), function(x) x@nugget.flag)
setMethod("nuggetflag", signature(x = "covAffineScaling"), function(x) x@nugget.flag)
setMethod("nuggetflag", signature(x = "covScaling"), function(x) x@nugget.flag)
setMethod("nuggetflag", signature(x = "covUser"), function(x) x@nugget.flag)

