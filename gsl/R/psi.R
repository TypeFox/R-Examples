"psi_int" <- function(n, give=FALSE, strict=TRUE){
  attr <- attributes(n)
  n.vec <- as.vector(n)
  jj <- .C("psi_int",
           as.integer(n.vec),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
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

"psi" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("psi",
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

"psi_1piy" <- function(y, give=FALSE, strict=TRUE){
  attr <- attributes(y)
  y.vec <- as.vector(y)
  jj <- .C("psi_1piy",
           as.double(y.vec),
           as.integer(length(y.vec)),
           val=as.double(y.vec),
           err=as.double(y.vec),
           status=as.integer(0*y.vec),
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

"psi_1_int" <- function(n, give=FALSE, strict=TRUE){
  attr <- attributes(n)
  n.vec <- as.vector(n)
  jj <- .C("psi_1_int",
           as.integer(n.vec),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
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

"psi_1" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("psi_1",
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

"psi_n" <- function(m, x, give=FALSE, strict=TRUE){
  jj <- process.args(m,x)
  m.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  jj <- .C("psi_n",
           as.integer(m.vec),
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
