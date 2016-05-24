"hyperg_0F1" <- function(c, x, give=FALSE, strict=TRUE){
  jj <- process.args(c,x)
  c.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("hyperg_0F1_e",
           as.double(c.vec),
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

"hyperg_1F1_int" <- function(m, n, x, give=FALSE, strict=TRUE){
  jj <- process.args(m,n,x)
  m.vec <- jj$arg1
  n.vec <- jj$arg2
  x.vec <- jj$arg3
  attr <- jj$attr

  jj <- .C("hyperg_1F1_int_e",
           as.integer(m.vec),
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



"hyperg_1F1" <- function(a, b, x, give=FALSE, strict=TRUE){
  jj <- process.args(a,b,x)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  x.vec <- jj$arg3
  attr <- jj$attr

  jj <- .C("hyperg_1F1_e",
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

"hyperg_U_int" <- function(m, n, x, give=FALSE, strict=TRUE){
  jj <- process.args(m,n,x)
  m.vec <- jj$arg1
  n.vec <- jj$arg2
  x.vec <- jj$arg3
  attr <- jj$attr

  jj <- .C("hyperg_U_int_e",
           as.integer(m.vec),
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


"hyperg_U" <- function(a, b, x, give=FALSE, strict=TRUE){
  jj <- process.args(a,b,x)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  x.vec <- jj$arg3
  attr <- jj$attr

  jj <- .C("hyperg_U_e",
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

"hyperg_2F1" <- function(a, b, c, x, give=FALSE, strict=TRUE){
  jj <- process.args(a, b, c, x)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  c.vec <- jj$arg3
  x.vec <- jj$arg4
  attr <- jj$attr

  jj <- .C("hyperg_2F1_e",
           as.double(a.vec),
           as.double(b.vec),
           as.double(c.vec),
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

"hyperg_2F1_conj" <- function(aR, aI, c, x, give=FALSE, strict=TRUE){
  jj <- process.args(aR,aI,c,x)
  aR.vec <- jj$arg1
  aI.vec <- jj$arg2
  c.vec <- jj$arg3
  x.vec <- jj$arg4
  attr <- jj$attr

  jj <- .C("hyperg_2F1_conj_e",
           as.double(aR.vec),
           as.double(aI.vec),
           as.double(c.vec),
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

"hyperg_2F1_renorm" <- function(a, b, c, x, give=FALSE, strict=TRUE){
  jj <- process.args(a,b,c,x)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  c.vec <- jj$arg3
  x.vec <- jj$arg4
  attr <- jj$attr

  jj <- .C("hyperg_2F1_renorm_e",
           as.double(a.vec),
           as.double(b.vec),
           as.double(c.vec),
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

"hyperg_2F1_conj_renorm" <- function(aR, aI, c, x, give=FALSE, strict=TRUE){
  jj <- process.args(aR,aI,c,x)
  aR.vec <- jj$arg1
  aI.vec <- jj$arg2
  c.vec <- jj$arg3
  x.vec <- jj$arg4
  attr <- jj$attr

  jj <- .C("hyperg_2F1_conj_renorm_e",
           as.double(aR.vec),
           as.double(aI.vec),
           as.double(c.vec),
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

"hyperg_2F0" <- function(a, b, x, give=FALSE, strict=TRUE){
  jj <- process.args(a,b,x)
  a.vec <- jj$arg1
  b.vec <- jj$arg2
  x.vec <- jj$arg3
  attr <- jj$attr

  jj <- .C("hyperg_2F0_e",
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

