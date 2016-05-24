identify.tri<-function(x,...)
  {
    if(!inherits(x,"tri"))
      stop("x must be of class \"tri\"")
    labels<-paste("(",x$x,",",x$y,")", sep ="")
    identify(x$x,x$y,labels=labels)
  }
