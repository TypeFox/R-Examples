##
## Methods for summarizing objects
## ===============================
##
## Method definition for objects of class "GoGARCH"
## Returns an object of class "Gosum" for which a show-method exists
##
setMethod(f = "summary", signature(object = "GoGARCH"), definition = function(object){
  name <- object@name
  method <- object@estby
  model <- object@garchf
  garchc <- lapply(object@models, function(x) x@fit$matcoef)
  ynames <- paste("y", 1:ncol(object@X), sep = "")
  names(garchc) <- paste("Component GARCH model of", ynames)
  garchc <- garchc
  Zinv <- solve(object@Z)
  gosum <- new("Gosum", name = name, method = method, model = model, garchc = garchc, Zinv = Zinv)
  return(gosum)
})
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "summary", signature(object = "Goestica"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "summary", signature(object = "Goestmm"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "summary", signature(object = "Goestnls"), definition = function(object){
  callNextMethod()
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "summary", signature(object = "Goestml"), definition = function(object){
  callNextMethod()
})

