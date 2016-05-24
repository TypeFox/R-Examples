setGeneric("contour", function(x, ...) standardGeneric("contour"))
setMethod(f = "contour", 
          signature = c(x = ".UD"), 
          definition = function(x, ...) {  
		  x<-getVolumeUD(x)
            callNextMethod()
          })


setMethod(f = "contour", 
          signature = c(x = ".UDStack"), 
          definition = function(x, ...){  
		  graphics::par(mfrow=rep(ceiling(sqrt(nlayers(x))), 2))
            lapply(split(x), contour, ...) 
          })
