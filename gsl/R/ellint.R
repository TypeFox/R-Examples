"ellint_Kcomp" <- function(k, mode=0, give=FALSE, strict=TRUE){
  attr <- attributes(k)
  if(length(mode)>1){stop("length of mode must be 1")}
  k.vec <- as.vector(k)
  jj <- .C("ellint_Kcomp_e",
           as.double(k.vec),
           as.integer(length(k.vec)),
           as.integer(mode),
           val=as.double(k.vec),
           err=as.double(k.vec),
           status=as.integer(0*k.vec),
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

"ellint_Ecomp" <- function(k, mode=0, give=FALSE, strict=TRUE){
  attr <- attributes(k)
  if(length(mode)>1){stop("length of mode must be 1")}
  k.vec <- as.vector(k)
  jj <- .C("ellint_Ecomp_e",
           as.double(k.vec),
           as.integer(length(k.vec)),
           as.integer(mode),
           val=as.double(k.vec),
           err=as.double(k.vec),
           status=as.integer(0*k.vec),
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

"ellint_F" <- function(phi, k, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}

  jj <- process.args(phi,k)
  phi.vec <- jj$arg1
  k.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("ellint_F_e",
           as.double(phi.vec),
           as.double(k.vec),
           as.integer(length(k.vec)),
           as.integer(mode),
           val=as.double(k.vec),
           err=as.double(k.vec),
           status=as.integer(0*k.vec),
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

 

"ellint_E" <- function(phi, k, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}

  jj <- process.args(phi,k)
  phi.vec <- jj$arg1
  k.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("ellint_E_e",
           as.double(phi.vec),
           as.double(k.vec),
           as.integer(length(k.vec)),
           as.integer(mode),
           val=as.double(k.vec),
           err=as.double(k.vec),
           status=as.integer(0*k.vec),
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

"ellint_P" <- function(phi, k, n, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}

  jj <- process.args(phi,k,n)
  phi.vec <- jj$arg1
  k.vec <- jj$arg2
  n.vec <- jj$arg3
  attr <- jj$attr
  
  jj <- .C("ellint_P_e",
           as.double(phi.vec),
           as.double(k.vec),
           as.double(n.vec),
           as.integer(length(k.vec)),
           as.integer(mode),
           val=as.double(k.vec),
           err=as.double(k.vec),
           status=as.integer(0*k.vec),
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

"ellint_D" <- function(phi, k, n, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}

  jj <- process.args(phi,k,n)
  phi.vec <- jj$arg1
  k.vec <- jj$arg2
  n.vec <- jj$arg3
  attr <- jj$attr
  
  jj <- .C("ellint_D_e",
           as.double(phi.vec),
           as.double(k.vec),
           as.double(n.vec),
           as.integer(length(k.vec)),
           as.integer(mode),
           val=as.double(k.vec),
           err=as.double(k.vec),
           status=as.integer(0*k.vec),
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

"ellint_RC" <- function(x, y, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}

  jj <- process.args(x,y)
  x.vec <- jj$arg1
  y.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("ellint_RC_e",
           as.double(x.vec),
           as.double(y.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(x.vec),
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

"ellint_RD" <- function(x, y, z, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}

  jj <- process.args(x,y,z)
  x.vec <- jj$arg1
  y.vec <- jj$arg2
  z.vec <- jj$arg3
  attr <- jj$attr
  
  jj <- .C("ellint_RD_e",
           as.double(x.vec),
           as.double(y.vec),
           as.double(z.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(x.vec),
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
"ellint_RF" <- function(x, y, z, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}

  jj <- process.args(x,y,z)
  x.vec <- jj$arg1
  y.vec <- jj$arg2
  z.vec <- jj$arg3
  attr <- jj$attr
  
  jj <- .C("ellint_RF_e",
           as.double(x.vec),
           as.double(y.vec),
           as.double(z.vec),
           as.integer(length(x.vec)),
           as.integer(mode),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(x.vec),
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


"ellint_RJ" <- function(x, y, z, p, mode=0, give=FALSE, strict=TRUE){
  if(length(mode)>1){stop("length of mode must be 1")}
  if(length(p)>1){stop("length of p must be 1")}

  jj <- process.args(x,y,z)
  x.vec <- jj$arg1
  y.vec <- jj$arg2
  z.vec <- jj$arg3
  attr <- jj$attr
  
  jj <- .C("ellint_RJ_e",
           as.double(x.vec),
           as.double(y.vec),
           as.double(z.vec),
           as.double(p),
           as.integer(length(x.vec)),
           as.integer(mode),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(x.vec),
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
