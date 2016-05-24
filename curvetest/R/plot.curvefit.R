plot.curvefit <-
structure(function(x, y=NULL, add=F, get.data=TRUE, ...){
      fits<-x
      y<-fits$data.model$y
      x<-fits$data.model$x
      if(!add) plot(fits$myx[o<-order(fits$myx)], fits$fitted[o],   ylim=range(y), xlab=fits$xlab, ylab=fits$ylab, 
      type="l",  ...) else 
      lines(fits$myx[o<-order(fits$myx)], fits$fitted[o], type="l",  ...) 
      #lines(x[o<-order(x)], fits$fits[o], col=3)
      if(get.data) points(x,  y, ...)
 }, modifiers = "public")
