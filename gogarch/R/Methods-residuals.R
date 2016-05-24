##
## Methods for returning the residuals of the GO-GARCH model
## =========================================================
##
## Method definition for objects of class "GoGARCH"
##
setMethod(f = "residuals", signature(object = "GoGARCH"), definition = function(object){
  m <- ncol(object@X)
  n <- nrow(object@X)
  svd <- lapply(object@H, svd)
  Vsqrinv <- lapply(svd, function(x) x$u %*% diag(1 / sqrt(x$d)) %*% t(x$u))
  resm <- matrix(0.0, nrow = n, ncol = m)
  for(i in 1:n){
    resm[i, ] <- Vsqrinv[[i]] %*% object@X[i, ]
  }
  colnames(resm) <- paste(colnames(object@X), "resid", sep = ".")
  rownames(resm) <- rownames(object@X)
  resm <- as.ts(resm)
  return(resm)
})
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "residuals", signature(object = "Goestica"), definition = function(object, standardize = FALSE){
  callNextMethod()
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "residuals", signature(object = "Goestmm"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "residuals", signature(object = "Goestnls"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "residuals", signature(object = "Goestml"), definition = function(object){
  callNextMethod()
})
