"EMPIRcopdf" <-
function(para=NULL, ...) {
  if(length(names(para)) != 2) {
    warning("a data.frame having only two columns is required")
    return(NULL)
  }

  uobs <- para[,1];
  vobs <- para[,2];
  n <- length(uobs);
  
  m <- length(uobs);
  if(m != length(vobs)) {
    warning("lengths of u and v are not identical")
    return(NULL) 
  }
  empcop <- EMPIRcop(uobs,vobs, para=para, ...)
  return(data.frame(u=uobs, v=vobs, empcop=empcop))
}


