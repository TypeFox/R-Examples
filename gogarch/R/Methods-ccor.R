##
## Methods for returning the conditional correlations
## ==================================================
##
## Method definition for objects of class "GoGARCH"
##
setMethod(f = "ccor", signature(object = "GoGARCH"), definition = function(object){
  m <- ncol(object@X)
  d <- m * (m - 1) / 2
  n <- nrow(object@X)
  cnames <- colnames(object@X)
  ccor <- matrix(c(unlist(lapply(object@H, function(x) cov2cor(x)[lower.tri(x)]))), ncol = d, nrow = n, byrow = TRUE)
  ngrid <- data.frame(expand.grid(cnames, cnames), stringsAsFactors = FALSE)
  mgrid <- paste(ngrid[, 1], ngrid[, 2], sep = " & ")
  mgrid <- matrix(mgrid, nrow = m, ncol = m)
  names <- mgrid[lower.tri(mgrid)]
  colnames(ccor) <- names
  rownames(ccor) <- rownames(object@X)
  ccor <- as.ts(ccor)
  return(ccor)
})
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "ccor", signature(object = "Goestica"), definition = function(object){
  ccor(as(object, "GoGARCH"))
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "ccor", signature(object = "Goestmm"), definition = function(object){
  ccor(as(object, "GoGARCH"))
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "ccor", signature(object = "Goestnls"), definition = function(object){
  ccor(as(object, "GoGARCH"))
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "ccor", signature(object = "Goestml"), definition = function(object){
  ccor(as(object, "GoGARCH"))
})
##
## Method definition for objects of class "Gopredict"
## "Gopredict" objects are returned by method "predict"
##
setMethod(f = "ccor", signature(object = "Gopredict"), definition = function(object){
  m <- ncol(object@Xf)
  d <- m * (m - 1) / 2
  n <- nrow(object@Xf)
  cnames <- colnames(object@Xf)
  ccor <- matrix(c(unlist(lapply(object@Hf, function(x) cov2cor(x)[lower.tri(x)]))), ncol = d, nrow = n, byrow = TRUE)
  ngrid <- data.frame(expand.grid(cnames, cnames), stringsAsFactors = FALSE)
  mgrid <- paste(ngrid[, 1], ngrid[, 2], sep = " & ")
  mgrid <- matrix(mgrid, nrow = m, ncol = m)
  names <- mgrid[lower.tri(mgrid)]
  colnames(ccor) <- names
  rownames(ccor) <- rownames(object@Xf)
  return(ccor)
})

