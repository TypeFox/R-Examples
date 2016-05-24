##
## Methods for obtaining predictions from GO-GARCH models
## ======================================================
##
## Method definition for objects of class "GoGARCH"
## The method "predict" returns an object of class "Gopredict" for
## which a "show" method exists.
##
setMethod(f = "predict", signature(object = "GoGARCH"), definition = function(object, n.ahead = 10, ...){
  n.ahead <- abs(as.integer(n.ahead))
  m <- ncol(object@X)
  n <- nrow(object@X)
  Z <- object@Z
  delta <- object@models[[1]]@fit$params$params["delta"]
  predictions <- lapply(object@models, predict, n.ahead = n.ahead)
  mean.pred.y <- matrix(unlist(lapply(predictions, function(x) x[, 1])), ncol = m)
  mean.pred.x <- mean.pred.y %*% Z
  rownames(mean.pred.x) <- seq(from = 1, to = n.ahead) + n
  colnames(mean.pred.x) <- paste(colnames(object@X), ".f", sep = "")
  h.pred.y <- matrix(unlist(lapply(predictions, function(x) x[, 3]^delta)), ncol = m)
  H.pred.y <- data.frame(t(h.pred.y))
  H.pred.x <- lapply(H.pred.y, function(x) Z %*% diag(x) %*% t(Z))
  names(H.pred.x) <- rownames(mean.pred.x)
  fcst <- new("Gopredict", Hf = H.pred.x, Xf = mean.pred.x, CGARCHF = predictions)
  return(fcst)
})
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "predict", signature(object = "Goestica"), definition = function(object, n.ahead = 10, ...){
  callNextMethod()
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "predict", signature(object = "Goestmm"), definition = function(object, n.ahead = 10, ...){
  callNextMethod()
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "predict", signature(object = "Goestnls"), definition = function(object, n.ahead = 10, ...){
  callNextMethod()
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "predict", signature(object = "Goestml"), definition = function(object, n.ahead = 10, ...){
  callNextMethod()
})
