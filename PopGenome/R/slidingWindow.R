setClass("slidingWindow", representation=representation(
      width		=    	"integer",
      jump		= 	"integer",
      type 		= 	"vector"),#, #erstes Vektorelement = Fenstertyp, zweite Vektorelement = Zusatzoption
      prototype 	= prototype( width = as.integer(50), jump = as.integer(50), type = c(1,0)
))

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol # R-Doku (?integer) statt is.integer(..)

setValidity("slidingWindow", function(object) {
      if(!is.wholenumber(object@width)) {
	    return("Window size must be integer!")
      }

      if(!is.wholenumber(object@jump)) {
	    return("Jump size must be integer!")
      }

      if(object@width < 1) {
	    return("Window size must be greater than 0!")
      }

      if(object@jump < 1) {
	    return("Jump size must be greater than 0!")
      }

      if(!(object@type[1] %in% c(1:4) )) {
  	    return("Window type must be between 1 and 4!")
      }

      if(object@type[2] < 0) {
  	    return("Window type option must be at least 0!")
      }

      # ansonsten:
      return(TRUE)
})