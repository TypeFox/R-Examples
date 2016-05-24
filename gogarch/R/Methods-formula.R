##
## Methods for returning the formula of the GARCH model
## ====================================================
##
## Method definition for objects of class "GoGARCH"
##
setMethod("formula", signature(x = "GoGARCH"), function(x, ...)
  x@garchf
)
##
## Method definition for objects of class "Goestica"
##
setMethod("formula", signature(x = "Goestica"), function(x, ...)
  x@garchf
)
##
## Method definition for objects of class "Goestmm"
##
setMethod("formula", signature(x = "Goestmm"), function(x, ...)
  x@garchf
)
##
## Method definition for objects of class "Goestnls"
##
setMethod("formula", signature(x = "Goestnls"), function(x, ...)
  x@garchf
)
##
## Method definition for objects of class "Goestml"
##
setMethod("formula", signature(x = "Goestml"), function(x, ...)
  x@garchf
)
