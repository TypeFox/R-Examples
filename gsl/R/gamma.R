"gsl_sf_gamma" <- function(x,give=FALSE,strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("gamma_e",
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

"lngamma" <- function(x,give=FALSE,strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("lngamma_e",
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

"lngamma_sgn" <- function(x, give=FALSE,strict=TRUE){
  jj <- process.args(x)
  x.vec <- jj$arg1
  attr <- jj$attr
  
  jj <- .C("lngamma_sgn_e",
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec),
           sgn=as.double(x.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  sgn <- jj$sgn
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr
  attributes(status) <- attr
  attributes(sgn) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status,sgn=sgn))
  } else {
    return(list(val=val,sgn=sgn))
  }
}

"gammastar" <- function(x,give=FALSE,strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("gammastar_e",
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

"gammainv" <- function(x,give=FALSE,strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("gammainv_e",
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

"lngamma_complex" <- function(zr, zi=NULL, r.and.i=TRUE, give=FALSE, strict=TRUE){
  attr <- attributes(zr)
  if(is.null(zi)){
    zi <- as.vector(Im(zr))
    zr <- as.vector(Re(zr))
  } else {
    zi <- as.vector(zi)
    zr <- as.vector(zr)
  }
  if(length(zr) !=length(zi)){stop("zr and zi must be of the same dimensions")}
  
  jj <- .C("lngamma_complex_e",
           as.double(zr),
           as.double(zi),
           as.integer(length(zr)),
           val_lnr=as.double(zr),
           val_arg=as.double(zr),
           err_lnr=as.double(zr),
           err_arg=as.double(zr),
           status=as.integer(0*zr),
           PACKAGE="gsl"
           )
  val_lnr <- jj$val_lnr
  val_arg <- jj$val_arg

  err_lnr <- jj$err_lnr
  err_arg <- jj$err_arg

  status <- jj$status
  attributes(status) <- attr
  
  if(r.and.i){
#    val <- exp(val_lnr)*cos(val_arg) + 1i*exp(val_lnr)*sin(val_arg)
#    err <- exp(xerr_lnr)*cos(err_arg) + 1i*exp(err_lnr)*sin(err_arg)
    val <- val_lnr + 1i*val_arg
    err <- err_lnr + 1i*err_arg
    attributes(val) <- attr
    attributes(err) <- attr

    if(strict){
      val <- strictify(val,status)
    }
    
    if(give){
      return(list(val=val, err=err, status=status))
    } else {
      return(val)
    }
  } else {
    attributes(val_lnr) <- attr
    attributes(val_arg) <- attr
    attributes(err_lnr) <- attr
    attributes(err_arg) <- attr

    if(strict){
      val_lnr <- strictify(val_lnr,status)
      val_arg <- strictify(val_arg,status)
    }

    if(give){
      return(list(val_lnr=val_lnr, val_arg=val_arg, err_lnr=err_lnr,err_arg=err_arg, status=status))
    } else {
      return(list(val_lnr=val_lnr, val_arg=val_arg))
    }
  } 
}

"taylorcoeff" <- function(n, x ,give=FALSE,strict=TRUE){
  jj <- process.args(n,x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  jj <- .C("taylorcoeff_e",
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
"fact" <- function(n, give=FALSE,strict=TRUE){
  n.vec <- as.vector(n)
  attr <- attributes(n)
  jj <- .C("fact_e",
           as.integer(n),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attributes(n)
  attributes(err) <- attributes(n)
  attributes(status) <- attributes(n)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"doublefact" <- function(n, give=FALSE,strict=TRUE){
  n.vec <- as.vector(n)
  attr <- attributes(n)
  jj <- .C("doublefact_e",
           as.integer(n),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attributes(n)
  attributes(err) <- attributes(n)
  attributes(status) <- attributes(n)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"lnfact" <- function(n, give=FALSE,strict=TRUE){
  n.vec <- as.vector(n)
  attr <- attributes(n)
  jj <- .C("lnfact_e",
           as.integer(n),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attributes(n)
  attributes(err) <- attributes(n)
  attributes(status) <- attributes(n)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"lndoublefact" <- function(n, give=FALSE,strict=TRUE){
  n.vec <- as.vector(n)
  attr <- attributes(n)
  jj <- .C("lndoublefact_e",
           as.integer(n),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attributes(n)
  attributes(err) <- attributes(n)
  attributes(status) <- attributes(n)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"gsl_sf_choose" <- function(n, m, give=FALSE,strict=TRUE){
  jj <- process.args(n,m)
  n.vec <- jj$arg1
  m.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("choose_e",
           as.integer(n.vec),
           as.integer(m.vec),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attributes(n)
  attributes(err) <- attributes(n)
  attributes(status) <- attributes(n)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}
"lnchoose" <- function(n, m, give=FALSE,strict=TRUE){
  jj <- process.args(n,m)
  n.vec <- jj$arg1
  m.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("lnchoose_e",
           as.integer(n.vec),
           as.integer(m.vec),
           as.integer(length(n.vec)),
           val=as.double(n.vec),
           err=as.double(n.vec),
           status=as.integer(0*n.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attributes(n)
  attributes(err) <- attributes(n)
  attributes(status) <- attributes(n)

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

"poch" <- function(a, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,x)
  a.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("poch_e",
           as.double(a.vec),
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

"lnpoch" <- function(a, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,x)
  a.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("lnpoch_e",
           as.double(a.vec),
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

"lnpoch_sgn" <- function(a, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,x)
  a.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("lnpoch_sgn_e",
           as.double(a.vec),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec),
           sgn=as.double(x.vec),
           PACKAGE="gsl"
           )
  val <- jj$val
  err <- jj$err
  sgn <- jj$sgn
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr
  attributes(status) <- attr
  attributes(sgn) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status,sgn=sgn))
  } else {
    return(list(val=val,sgn=sgn))
  }
}

"pochrel" <- function(a, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,x)
  a.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("pochrel_e",
           as.double(a.vec),
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

"gamma_inc_Q" <- function(a, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,x)
  a.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("gamma_inc_Q_e",
           as.double(a.vec),
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

"gamma_inc_P" <- function(a, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,x)
  a.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("gamma_inc_P_e",
           as.double(a.vec),
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

"gamma_inc" <- function(a, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,x)
  a.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("gamma_inc_e",
           as.double(a.vec),
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

"gsl_sf_beta" <- function(a, b, give=FALSE,strict=TRUE){
  jj <- process.args(a,b)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("beta_e",
           as.double(a.vec),
           as.double(b.vec),
           as.integer(length(b.vec)),
           val=as.double(b.vec),
           err=as.double(b.vec),
           status=as.integer(0*b.vec),
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

"lnbeta" <- function(a, b, give=FALSE,strict=TRUE){
  jj <- process.args(a,b)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("lnbeta_e",
           as.double(a.vec),
           as.double(b.vec),
           as.integer(length(b.vec)),
           val=as.double(b.vec),
           err=as.double(b.vec),
           status=as.integer(0*b.vec),
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

"beta_inc" <- function(a, b, x, give=FALSE,strict=TRUE){
  jj <- process.args(a,b,x)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  x.vec <- jj$arg3
  attr <- jj$attr
  
  jj <- .C("beta_inc_e",
           as.double(a.vec),
           as.double(b.vec),
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
