"coupling_3j" <- function(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc, give=FALSE, strict=TRUE){
  jj <- process.args(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc)
  ja.vec <- jj$arg1
  jb.vec <- jj$arg2
  jc.vec <- jj$arg3
  ma.vec <- jj$arg4
  mb.vec <- jj$arg5
  mc.vec <- jj$arg6
  attr <- jj$attr

  jj <- .C("coupling_3j",
           as.integer(ja.vec), 
           as.integer(jb.vec), 
           as.integer(jc.vec), 
           as.integer(ma.vec), 
           as.integer(mb.vec), 
           as.integer(mc.vec), 
           as.integer(length(ja.vec)),
           val=as.double(ja.vec),
           err=as.double(ja.vec),
           status=as.integer(0*ja.vec),
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

"coupling_6j" <- function(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, give=FALSE, strict=TRUE){
  jj <- process.args(two_ja, two_jb, two_jc, two_jd, two_je, two_jf)
  ja.vec <- jj$arg1
  jb.vec <- jj$arg2
  jc.vec <- jj$arg3
  jd.vec <- jj$arg4
  je.vec <- jj$arg5
  jf.vec <- jj$arg6
  attr <- jj$attr

  jj <- .C("coupling_6j",
           as.integer(ja.vec), 
           as.integer(jb.vec), 
           as.integer(jc.vec), 
           as.integer(jd.vec), 
           as.integer(je.vec), 
           as.integer(jf.vec), 
           as.integer(length(ja.vec)),
           val=as.double(ja.vec),
           err=as.double(ja.vec),
           status=as.integer(0*ja.vec),
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

"coupling_9j" <- function(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji, give=FALSE, strict=TRUE){
  jj <- process.args(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji)
  ja.vec <- jj$arg1
  jb.vec <- jj$arg2
  jc.vec <- jj$arg3
  jd.vec <- jj$arg4
  je.vec <- jj$arg5
  jf.vec <- jj$arg6
  jg.vec <- jj$arg7
  jh.vec <- jj$arg8
  ji.vec <- jj$arg9
  attr <- jj$attr

  jj <- .C("coupling_9j",
           as.integer(ja.vec), 
           as.integer(jb.vec), 
           as.integer(jc.vec), 
           as.integer(jd.vec), 
           as.integer(je.vec), 
           as.integer(jf.vec), 
           as.integer(jg.vec), 
           as.integer(jh.vec), 
           as.integer(ji.vec), 
           as.integer(length(ja.vec)),
           val=as.double(ja.vec),
           err=as.double(ja.vec),
           status=as.integer(0*ja.vec),
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

