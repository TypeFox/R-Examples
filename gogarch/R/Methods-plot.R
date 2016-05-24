##
## Methods for plotting GO-GARCH models
## ====================================
##
## Method definition for objects of class "GoGARCH"
##
setMethod(f = "plot", signature(x = "GoGARCH", y = "missing"), definition = function(x, main = NULL, ...){
  if(is.null(main)){
    main <- paste("Conditional correlations of", x@name)
  }
  x <- ccor(x)
  plot(x, main = main, ...)
})
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "plot", signature(x = "Goestica", y = "missing"), definition = function(x, main = NULL, ...){
  x <- as(x, "GoGARCH")
  plot(x, main = main, ...)
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "plot", signature(x = "Goestmm", y = "missing"), definition = function(x, main = NULL, ...){
  x <- as(x, "GoGARCH")
  plot(x, main = main, ...)
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "plot", signature(x = "Goestnls", y = "missing"), definition = function(x, main = NULL, ...){
  x <- as(x, "GoGARCH")
  plot(x, main = main, ...)
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "plot", signature(x = "Goestml", y = "missing"), definition = function(x, main = NULL, ...){
  x <- as(x, "GoGARCH")
  plot(x, main = main, ...)
})
