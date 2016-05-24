"expint_E1" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("expint_E1_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"expint_E2" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("expint_E2_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"expint_En" <- function(n, x, give=FALSE, strict=TRUE){
  jj <- process.args(n, x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("expint_En_e",
           as.integer(n.vec),
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

"expint_Ei" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("expint_Ei_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}
"Shi" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("Shi_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"Chi" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("Chi_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"expint_3" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("expint_3_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"Si" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("Si_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"Ci" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("Ci_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"atanint" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  jj <- .C("atanint_e",
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
  attributes(val) <- attributes(x)
  attributes(err) <- attributes(x)
  attributes(status) <- attributes(x)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}
