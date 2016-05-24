"bessel_J0" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_J0_e",
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
   
"bessel_J1" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_J1_e",
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

"bessel_Jn" <- function(n, x, give=FALSE, strict=TRUE){
  jj <- process.args(n, x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_Jn_e",
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

"bessel_Jn_array" <- function(nmin,nmax, x, give=FALSE,strict=TRUE){
  if(length(nmin)>1){stop("nmin should be of length 1")}
  if(length(nmax)>1){stop("nmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(nmax-nmin+1))
  jj <- .C("bessel_Jn_array_e",
           as.integer(nmin),
           as.integer(nmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(nmax-nmin+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_Y0" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_Y0_e",
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
           
"bessel_Y1" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_Y1_e",
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

"bessel_Yn" <- function(n, x, give=FALSE, strict=TRUE){
  jj <- process.args(n,x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  jj <- .C("bessel_Yn_e",
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

"bessel_Yn_array" <- function(nmin,nmax, x, give=FALSE,strict=TRUE){
  if(length(nmin)>1){stop("nmin should be of length 1")}
  if(length(nmax)>1){stop("nmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(nmax-nmin+1))
  jj <- .C("bessel_Yn_array_e",
           as.integer(nmin),
           as.integer(nmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(nmax-nmin+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_I0" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_I0_e",
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
   
"bessel_I1" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_I1_e",
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

"bessel_In" <- function(n, x, give=FALSE, strict=TRUE){
  jj <- process.args(n, x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_In_e",
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

"bessel_In_array" <- function(nmin,nmax, x, give=FALSE,strict=TRUE){
  if(length(nmin)>1){stop("nmin should be of length 1")}
  if(length(nmax)>1){stop("nmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(nmax-nmin+1))
  jj <- .C("bessel_In_array_e",
           as.integer(nmin),
           as.integer(nmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(nmax-nmin+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_I0_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_I0_scaled_e",
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
   
"bessel_I1_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_I1_scaled_e",
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

"bessel_In_scaled" <- function(n, x, give=FALSE, strict=TRUE){
  jj <- process.args(n,x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("bessel_In_scaled_e",
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

"bessel_In_scaled_array" <- function(nmin,nmax, x, give=FALSE,strict=TRUE){
  if(length(nmin)>1){stop("nmin should be of length 1")}
  if(length(nmax)>1){stop("nmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(nmax-nmin+1))
  jj <- .C("bessel_In_scaled_array_e",
           as.integer(nmin),
           as.integer(nmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(nmax-nmin+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_K0" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_K0_e",
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
   
"bessel_K1" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_K1_e",
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

"bessel_Kn" <- function(n, x, give=FALSE, strict=TRUE){
  jj <- process.args(n, x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("bessel_Kn_e",
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

"bessel_Kn_array" <- function(nmin,nmax, x, give=FALSE,strict=TRUE){
  if(length(nmin)>1){stop("nmin should be of length 1")}
  if(length(nmax)>1){stop("nmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(nmax-nmin+1))
  jj <- .C("bessel_Kn_array_e",
           as.integer(nmin),
           as.integer(nmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(nmax-nmin+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_K0_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_K0_scaled_e",
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
   
"bessel_K1_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_K1_scaled_e",
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

"bessel_Kn_scaled" <- function(n, x, give=FALSE, strict=TRUE){
  jj <- process.args(n, x)
  n.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_Kn_scaled_e",
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

"bessel_Kn_scaled_array" <- function(nmin,nmax, x, give=FALSE,strict=TRUE){
  if(length(nmin)>1){stop("nmin should be of length 1")}
  if(length(nmax)>1){stop("nmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  nmin <- as.integer(nmin)
  nmax <- as.integer(nmax)
  x.out <- rep(x.vec,(nmax-nmin+1))
  jj <- .C("bessel_Kn_scaled_array_e",
           as.integer(nmin),
           as.integer(nmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(nmax-nmin+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_j0" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_j0_e",
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
   
"bessel_j1" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_j1_e",
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

"bessel_j2" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_j2_e",
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

"bessel_jl" <- function(l, x, give=FALSE, strict=TRUE){
  jj <- process.args(l, x)
  l.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr


  jj <- .C("bessel_jl_e",
           as.integer(l.vec),
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

"bessel_jl_array" <- function(lmax, x, give=FALSE,strict=TRUE){
  if(length(lmax)>1){stop("lmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(lmax+1))
  jj <- .C("bessel_jl_array_e",
           as.integer(lmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(lmax+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_jl_steed_array" <- function(lmax, x, give=FALSE,strict=TRUE){
  if(length(lmax)>1){stop("lmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(lmax+1))
  jj <- .C("bessel_jl_steed_array_e",
           as.integer(lmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(lmax+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_y0" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_y0_e",
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
   
"bessel_y1" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_y1_e",
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

"bessel_y2" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_y2_e",
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

"bessel_yl" <- function(l, x, give=FALSE, strict=TRUE){
  jj <- process.args(l, x)
  l.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$arg3
  
  jj <- .C("bessel_yl_e",
           as.integer(l.vec),
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

"bessel_yl_array" <- function(lmax, x, give=FALSE, strict=TRUE){
  if(length(lmax)>1){stop("lmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(lmax+1))
  jj <- .C("bessel_yl_array_e",
           as.integer(lmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(lmax+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_i0_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_i0_scaled_e",
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
   
"bessel_i1_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_i1_scaled_e",
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

"bessel_i2_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_i2_scaled_e",
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

"bessel_il_scaled" <- function(l, x, give=FALSE, strict=TRUE){
  jj <- process.args(l, x)
  l.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("bessel_il_scaled_e",
           as.integer(l.vec),
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

"bessel_il_scaled_array" <- function(lmax, x, give=FALSE,strict=TRUE){
  if(length(lmax)>1){stop("lmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec,(lmax+1))
  jj <- .C("bessel_il_scaled_array_e",
           as.integer(lmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(lmax+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_k0_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_k0_scaled_e",
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
   
"bessel_k1_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_k1_scaled_e",
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

"bessel_k2_scaled" <- function(x, give=FALSE, strict=TRUE){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("bessel_k2_scaled_e",
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

"bessel_kl_scaled" <- function(l, x, give=FALSE, strict=TRUE){
  jj <- process.args(l, x)
  l.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_kl_scaled_e",
           as.integer(l.vec),
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

"bessel_kl_scaled_array" <- function(lmax, x, give=FALSE,strict=TRUE){
  if(length(lmax)>1){stop("lmax should be of length 1")}
  x.vec <- as.vector(x)
  attr <- attributes(x)
  x.out <- rep(x.vec, (lmax+1))
  jj <- .C("bessel_kl_scaled_array_e",
           as.integer(lmax),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(lmax+1, length(x.vec))
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,status=status))
  } else {
    return(val)
  }  
}  

"bessel_Jnu" <- function (nu, x, give = FALSE, strict = TRUE){
  jj <- process.args(nu, x)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_Jnu_e",
           as.double(nu.vec),
           as.double(x.vec), 
           as.integer(length(x.vec)),
           val = as.double(x.vec),
           err = as.double(x.vec), 
           status = as.integer(0*x.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr
  
  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_sequence_Jnu" <- function (nu, v, mode = 0, give = FALSE, strict=TRUE){
  if (length(nu) > 1 ) {
    stop("nu should be of length 1")
  }
  if (length(mode) > 1 ) {
    stop("mode should be of length 1")
  }
  v.vec <- as.vector(v)
  if(any(v.vec<0)){stop("all elements of v must be positive")}
  if( !all(diff(v.vec)>0) ){stop("elements of v must be sorted in increasing order")}
     
  jj <- .C("bessel_sequence_Jnu_e",
           as.double(nu),
           val=as.double(v.vec), 
           as.integer(length(v.vec)),
           as.integer(mode),
           status=as.integer(nu),
           PACKAGE = "gsl"
           )
  val <- jj$val
  status <- jj$status
  attributes(val) <- attributes(v)
  attributes(status) <- attributes(v)

  if(strict){
    if(status>0){val[] <- NA}
  }
  
  if(give){
    return(list(val=val, status=status))
  } else {    
    return(val)
  }
}

"bessel_Ynu" <- function (nu, x, give = FALSE, strict = TRUE){
  jj <- process.args(nu, x)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_Ynu_e",
           as.double(nu),
           as.double(x.vec), 
           as.integer(length(x.vec)),
           val = as.double(x.vec),
           err = as.double(x.vec), 
           status = as.integer(0*x.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  }
  else {
    return(val)
  }
}

"bessel_Inu" <- function (nu, x, give = FALSE, strict = TRUE){
  jj <- process.args(nu, x)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("bessel_Inu_e",
           as.double(nu.vec),
           as.double(x.vec), 
           as.integer(length(x.vec)),
           val = as.double(x.vec),
           err = as.double(x.vec), 
           status = as.integer(0*x.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_Inu_scaled" <- function (nu, x, give = FALSE, strict = TRUE){
  jj <- process.args(nu, x)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_Inu_scaled_e",
           as.double(nu.vec),
           as.double(x.vec), 
           as.integer(length(x.vec)),
           val = as.double(x.vec),
           err = as.double(x.vec), 
           status = as.integer(0*x.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_Knu" <- function (nu, x, give = FALSE, strict = TRUE){
  jj <- process.args(nu, x)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_Knu_e",
           as.double(nu.vec),
           as.double(x.vec), 
           as.integer(length(x.vec)),
           val = as.double(x.vec),
           err = as.double(x.vec), 
           status = as.integer(0*x.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_lnKnu" <- function (nu, x, give = FALSE, strict = TRUE){
  jj <- process.args(nu, x)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_lnKnu_e",
           as.double(nu.vec),
           as.double(x.vec), 
           as.integer(length(x.vec)),
           val = as.double(x.vec),
           err = as.double(x.vec), 
           status = as.integer(0*x.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_Knu_scaled" <- function (nu, x, give = FALSE, strict = TRUE){
  jj <- process.args(nu, x)
  nu.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("bessel_Knu_scaled_e",
           as.double(nu.vec),
           as.double(x.vec), 
           as.integer(length(x.vec)),
           val = as.double(x.vec),
           err = as.double(x.vec), 
           status = as.integer(0*x.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_zero_J0" <- function (s, give = FALSE, strict = TRUE){
  s.vec <- as.vector(s)
  attr <- attributes(s)
  jj <- .C("bessel_zero_J0_e",
           as.integer(s.vec), 
           as.integer(length(s.vec)),
           val = as.double(s.vec),
           err = as.double(s.vec), 
           status = as.integer(s.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_zero_J1" <- function (s, give = FALSE, strict = TRUE){
  s.vec <- as.vector(s)
  attr <- attributes(s)
  jj <- .C("bessel_zero_J1_e", as.integer(s.vec), 
           as.integer(length(s.vec)),
           val = as.double(s.vec),
           err = as.double(s.vec), 
           status = as.integer(s.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
    }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}

"bessel_zero_Jnu" <- function (nu, s, give = FALSE, strict = TRUE){
  jj <- process.args(nu,s)
  nu.vec <- jj$arg1
  s.vec <- jj$arg2
  attr <- jj$attr
  jj <- .C("bessel_zero_Jnu_e",
           as.double(nu.vec),
           as.integer(s.vec), 
           as.integer(length(s.vec)),
           val = as.double(s.vec),
           err = as.double(s.vec), 
           status = as.integer(s.vec),
           PACKAGE = "gsl"
           )
  val <- jj$val
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr  
  attributes(status) <- attr

  if (strict) {
    val <- strictify(val, status)
  }
  if (give) {
    return(list(val = val, err = err, status = status))
  } else {
    return(val)
  }
}
