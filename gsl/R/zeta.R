"zeta_int" <- function(n, give=FALSE, strict=TRUE){
  attr <- attributes(n)
  n.vec <- as.vector(n)
  jj <- .C("zeta_int",
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

"zeta" <- function(s, give=FALSE, strict=TRUE){
  attr <- attributes(s)
  s.vec <- as.vector(s)
  jj <- .C("zeta",
           as.double(s.vec),
           as.integer(length(s.vec)),
           val=as.double(s.vec),
           err=as.double(s.vec),
           status=as.integer(0*s.vec),
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

"zetam1_int" <- function(n, give=FALSE, strict=TRUE){
  attr <- attributes(n)
  n.vec <- as.vector(n)
  jj <- .C("zetam1_int",
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

"zetam1" <- function(s, give=FALSE, strict=TRUE){
  attr <- attributes(s)
  s.vec <- as.vector(s)
  jj <- .C("zetam1",
           as.double(s.vec),
           as.integer(length(s.vec)),
           val=as.double(s.vec),
           err=as.double(s.vec),
           status=as.integer(0*s.vec),
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

"hzeta" <- function(s, q, give=FALSE, strict=TRUE){
  jj <- process.args(s,q)
  s.vec <- jj$arg1
  q.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("hzeta",
           as.double(s.vec),
           as.double(q.vec),
           as.integer(length(s.vec)),
           val=as.double(s.vec),
           err=as.double(s.vec),
           status=as.integer(0*s.vec),
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

"eta_int" <- function(n, give=FALSE, strict=TRUE){
  attr <- attributes(n)
  n.vec <- as.vector(n)
  jj <- .C("eta_int",
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

"eta" <- function(s, give=FALSE, strict=TRUE){
  attr <- attributes(s)
  s.vec <- as.vector(s)
  jj <- .C("eta",
           as.double(s.vec),
           as.integer(length(s.vec)),
           val=as.double(s.vec),
           err=as.double(s.vec),
           status=as.integer(0*s.vec),
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

