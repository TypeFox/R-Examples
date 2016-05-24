"gsl_sf_log" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("log_e",
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

"log_abs" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("log_abs_e",
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

"complex_log" <- function(zr, zi=NULL, r.and.i=TRUE, give=FALSE, strict=TRUE){
  attr <- attributes(zr)
  if(is.null(zi)){
    zi <- as.vector(Im(zr))
    zr <- as.vector(Re(zr))
  } else {
    zi <- as.vector(zi)
    zr <- as.vector(zr)
  }
  if(length(zr) !=length(zi)){stop("zr and zi must be of the same dimensions")}
  
  jj <- .C("complex_log_e",
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

"log_1plusx" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("log_1plusx_e",
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

"log_1plusx_mx" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("log_1plusx_mx_e",
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
