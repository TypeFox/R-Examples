allele.names<-function(x)
  {
    retval  <- attr(x,"allele.names")
    if(is.null(retval))
      retval  <- x$allele.names
    return(retval)
  }




