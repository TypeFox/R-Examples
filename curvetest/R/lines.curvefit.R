lines.curvefit <-
structure(function(x, ...){
      fits<-x
      x<-fits$myx
      y<-fits$fitted
      lines(x[o<-order(x)], y[o],   ...)        
 }, modifiers = "public")
