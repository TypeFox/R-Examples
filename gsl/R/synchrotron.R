"synchrotron_1" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("synchrotron_1",
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr


  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
      return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}  

"synchrotron_2" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("synchrotron_2",
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
      return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}  

