"dilog" <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("dilog_e",
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

"complex_dilog" <- function(r, theta=NULL, give=FALSE, strict=TRUE){  
  if(is.null(theta)){
    attr <- attributes(r)
    r.vec <- as.vector(Mod(r))
    theta.vec <- as.vector(Arg(r))
  } else {
    jj <- process.args(r,theta)
    r.vec <- jj$arg1
    theta.vec <- jj$arg2
    attr <- jj$attr
  }
  
  jj <- .C("complex_dilog_e",
           as.double(r.vec),
           as.double(theta),
           as.integer(length(r)),
           val_re=as.double(r.vec),
           val_im=as.double(r.vec),
           err_re=as.double(r.vec),
           err_im=as.double(r.vec),
           status=as.integer(0*r.vec),
           PACKAGE="gsl"
           )
  
  val <- jj$val_re + 1i*jj$val_im
  err <- jj$err_re + 1i*jj$err_im
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr
  attributes(status) <- attr
  
  if(give){
    return(list(val=val, err=err, status=status))
  } else {
    return(val)
  }
}
