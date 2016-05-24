"erf" <- function(x, mode=0, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("erf_e",
           as.double(x.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
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

"erfc" <- function(x, mode=0, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("erfc_e",
           as.double(x.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
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

"log_erfc" <- function(x, mode=0, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("log_erfc_e",
           as.double(x.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
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

"log_erf_Z" <- function(x, mode=0, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("erf_e",
           as.double(x.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
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

"erf_Q" <- function(x, mode=0, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("erf_Q_e",
           as.double(x.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
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

"hazard" <- function(x, mode=0, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("hazard_e",
           as.double(x.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
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
