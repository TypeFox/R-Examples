"gegenpoly_1" <- function(lambda, x, give=FALSE, strict=TRUE){

  jj <- process.args(lambda,x)
  lambda.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("gegenpoly_1_e",
           as.double(lambda.vec),
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

"gegenpoly_2" <- function(lambda, x, give=FALSE, strict=TRUE){

  jj <- process.args(lambda,x)
  lambda.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("gegenpoly_2_e",
           as.double(lambda.vec),
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

"gegenpoly_3" <- function(lambda, x, give=FALSE, strict=TRUE){

  jj <- process.args(lambda,x)
  lambda.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("gegenpoly_3_e",
           as.double(lambda.vec),
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


"gegenpoly_n" <- function(n, lambda, x, give=FALSE, strict=TRUE){

  if(length(n)>1){stop("length of n should be 1")}
  jj <- process.args(lambda,x)
  lambda.single <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  jj <- .C("gegenpoly_n_e",
           as.integer(n),
           as.double(lambda.single),
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

"gegenpoly_array" <- function(nmax, lambda, x, give=FALSE,strict=TRUE){
  if(length(nmax)>1){stop("nmax should be of length 1")}
  jj <- process.args(lambda,x)
  lambda.single <- jj$arg1
  x.vec<- jj$arg2
  attr <- jj$attr
  x.out <- rep(x.vec,(nmax+1))
  
  jj <- .C("gegenpoly_array",
           as.integer(nmax),
           as.double(lambda.single),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(nmax+1 , length(x.vec))
  status <- jj$status
  attributes(status) <- attr
  err <- jj$err
  attributes(err) <- attr
  
  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }  
}  

