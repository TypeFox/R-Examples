"prepanel.dice" <-
function(x,y){
  xx <- ceiling(sqrt(length(x)))
  yy <- ceiling( length(x)/xx )
  return(list(ylim=c(-0.1,yy+0.1),xlim=c(-0.1,xx+0.1)) )
}

