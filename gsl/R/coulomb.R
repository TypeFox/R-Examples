"hydrogenicR_1" <- function(Z, r, give=FALSE, strict=TRUE){
  jj <- process.args(Z,r)
  Z.vec <- jj$arg1
  r.vec <- jj$arg2
  attr <- jj$attr
  jj <- .C("hydrogenicR_1",
           as.double(Z.vec),
           as.double(r.vec),
           as.integer(length(Z.vec)),
           val=as.double(Z.vec),
           err=as.double(Z.vec),
           status=as.integer(0*Z.vec),
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

"hydrogenicR" <- function(n, l, Z, r, give=FALSE, strict=TRUE){
  jj <- process.args(n,l,Z,r)
  n.vec <- jj$arg1
  l.vec <- jj$arg2
  Z.vec <- jj$arg3
  r.vec <- jj$arg4
  attr <- jj$attr
  jj <- .C("hydrogenicR",
           as.integer(n.vec),
           as.integer(l.vec),
           as.double(Z.vec),
           as.double(r.vec),
           as.integer(length(r.vec)),
           val=as.double(r.vec),
           err=as.double(r.vec),
           status=as.integer(0*r.vec),
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

"coulomb_wave_FG" <- function(eta, x, L_F, k, give=FALSE, strict=TRUE){
  jj <- process.args(eta, x, L_F, k)
  eta.vec <- jj$arg1
  x.vec <- jj$arg2
  L_F.vec <- jj$arg3
  k.vec <- jj$arg4
  attr <- jj$attr
  jj <- .C("coulomb_wave_FG",
           as.double(eta.vec),
           as.double(x.vec),
           as.double(L_F.vec),
           as.integer(k.vec),
           as.integer(length(eta.vec)),
           val_F=as.double(0*eta.vec),
           err_F=as.double(0*eta.vec),
           val_Fp=as.double(0*eta.vec),
           err_Fp=as.double(0*eta.vec),
           val_G=as.double(0*eta.vec),
           err_G=as.double(0*eta.vec),
           val_Gp=as.double(0*eta.vec),
           err_Gp=as.double(0*eta.vec),
           exp_F=as.double(0*eta.vec),
           exp_G=as.double(0*eta.vec),
           status=as.integer(0*eta.vec),
           PACKAGE="gsl"
           )
  val_F <- jj$val_F
  val_Fp <- jj$val_Fp
  val_G <- jj$val_G
  val_Gp <- jj$val_Gp

  err_F <- jj$err_F
  err_Fp <- jj$err_Fp
  err_G <- jj$err_Gp
  err_Gp <- jj$err_Gp

  status <- jj$status

  exp_F <- jj$exp_F
  exp_G <- jj$exp_G
  
  attributes(val_F) <- attr
  attributes(val_Fp) <- attr
  attributes(val_G) <- attr
  attributes(val_Gp) <- attr

  attributes(err_F) <- attr
  attributes(err_Fp) <- attr
  attributes(err_G) <- attr
  attributes(err_Gp) <- attr

  attributes(exp_F) <- attr
  attributes(exp_G) <- attr
  
  attributes(status) <- attr
  
  if(strict){
    val_F <- strictify(val_F,status)
    val_Fp <- strictify(val_Fp,status)
    val_G <- strictify(val_G,status)
    val_Gp <- strictify(val_Gp,status)

    err_F <- strictify(err_F,status)
    err_Fp <- strictify(err_Fp,status)
    err_G <- strictify(err_G,status)
    err_Gp <- strictify(err_Gp,status)

    exp_F <- strictify(exp_F,status)
    exp_G <- strictify(exp_G,status)
  }
  
  if(give){
      return(list(val_F=val_F,
                  val_Fp=val_Fp,
                  val_G=val_G,
                  val_Gp=val_Gp,
                  err_F=err_F,
                  err_Fp=err_Fp,
                  err_G=err_G,
                  err_Gp=err_Gp,
                  exp_F=exp_F,
                  exp_G=exp_G,
                  status=status
                  )
             )
      
  } else {
      return(list(val_F=val_F,
                  val_Fp=val_Fp,
                  val_G=val_G,
                  val_Gp=val_Gp,
                  exp_F=exp_F,
                  exp_G=exp_G
                  )
             )
  }
}  

"coulomb_wave_F_array" <- function(L_min, kmax, eta, x, give=FALSE,strict=TRUE){
  if(length(L_min)>1){stop("L_min should be of length 1")}
  if(length(kmax)>1){stop("kmax should be of length 1")}
  jj <- process.args(eta,x)
  eta.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  x.out <- rep(x.vec,(kmax+1))
  jj <- .C("coulomb_wave_F_array",
           as.double(L_min),
           as.integer(kmax),
           as.double(eta.vec),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           F_exp=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(kmax+1, length(x.vec))
  rownames(val) <- L_min:(L_min+kmax)
  colnames(val) <- names(x)
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val, F_exp=jj$F_exp, status=status))
  } else {
    return(val)
  }  
}  

"coulomb_wave_FG_array" <- function(L_min, kmax, eta, x, give=FALSE,strict=TRUE){
  if(length(L_min)>1){stop("L_min should be of length 1")}
  if(length(kmax)>1){stop("kmax should be of length 1")}
  jj <- process.args(eta,x)
  eta.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  x.out <- rep(x.vec,(kmax+1))
  jj <- .C("coulomb_wave_FG_array",
           as.double(L_min),
           as.integer(kmax),
           as.double(eta.vec),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val_F=as.double(x.out),
           val_G=as.double(x.out),
           F_exp=as.double(x.vec),
           G_exp=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val_F <- jj$val_F
  val_G <- jj$val_G
  F_exp <- jj$F_exp
  G_exp <- jj$G_exp
  dim(val_F) <- c(kmax+1, length(x.vec))
  dim(val_G) <- c(kmax+1, length(x.vec))

  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val_F <- strictify(val_F,status)
    val_G <- strictify(val_G,status)
  }
  
  if(give){
    return(list(val_F=val_F, val_G=val_G, F_exp=F_exp, status=status))
  } else {
    return(list(val_F=val_F, val_G=val_G))
  }  
}

"coulomb_wave_FGp_array" <- function(L_min, kmax, eta, x, give=FALSE,strict=TRUE){
  if(length(L_min)>1){stop("L_min should be of length 1")}
  if(length(kmax)>1){stop("kmax should be of length 1")}
  jj <- process.args(eta,x)
  eta.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  x.out <- rep(x.vec,(kmax+1))
  jj <- .C("coulomb_wave_FGp_array",
           as.double(L_min),
           as.integer(kmax),
           as.double(eta.vec),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val_F=as.double(x.out),
           val_Fp=as.double(x.out),
           val_G=as.double(x.out),
           val_Gp=as.double(x.out),
           F_exp=as.double(x.vec),
           G_exp=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val_F <- jj$val_F
  val_Fp <- jj$val_Fp
  val_G <- jj$val_G
  val_Gp <- jj$val_Gp
  F_exp <- jj$F_exp
  G_exp <- jj$G_exp
  dim(val_F) <- c(kmax+1, length(x.vec))
  dim(val_Fp) <- c(kmax+1, length(x.vec))
  dim(val_G) <- c(kmax+1, length(x.vec))
  dim(val_Gp) <- c(kmax+1, length(x.vec))

  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val_F <- strictify(val_F,status)
    val_Fp <- strictify(val_Fp,status)
    val_G <- strictify(val_G,status)
    val_Gp <- strictify(val_Gp,status)
  }
  
  if(give){
    return(list(val_F=val_F, val_Fp=val_Fp, val_G=val_G, val_Gp=val_Gp, F_exp=F_exp, status=status))
  } else {
    return(list(val_F=val_F, val_Fp=val_Fp, val_G=val_G, val_Gp=val_Gp))
  }  
}

"coulomb_wave_sphF_array" <- function(L_min, kmax, eta, x, give=FALSE,strict=TRUE){
  if(length(L_min)>1){stop("L_min should be of length 1")}
  if(length(kmax)>1){stop("kmax should be of length 1")}
  jj <- process.args(eta,x)
  eta.vec <- jj$arg1
  x.vec <- jj$arg2
  attr <- jj$attr
  
  x.out <- rep(x.vec,(kmax+1))
  jj <- .C("coulomb_wave_sphF_array",
           as.double(L_min),
           as.integer(kmax),
           as.double(eta.vec),
           as.double(x.vec),
           as.integer(length(x.vec)),
           val=as.double(x.out),
           F_exp=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(kmax+1, length(x.vec))
  rownames(val) <- L_min:(L_min+kmax)
  colnames(val) <- names(x)
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val, F_exp=jj$F_exp, status=status))
  } else {
    return(val)
  }  
}

"coulomb_CL" <- function(L, eta, give=FALSE, strict=TRUE){
  jj <- process.args(L,eta)
  L.vec <- jj$arg1
  eta.vec <- jj$arg2
  attr <- jj$attr
  jj <- .C("coulomb_CL",
           as.double(L.vec),
           as.double(eta.vec),
           as.integer(length(L.vec)),
           val=as.double(L.vec),
           err=as.double(L.vec),
           status=as.integer(0*L.vec),
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

"coulomb_CL_array" <- function(L_min, kmax, eta, give=FALSE,strict=TRUE){
  if(length(L_min)>1){stop("L_min should be of length 1")}
  if(length(kmax)>1){stop("kmax should be of length 1")}
  jj <- process.args(eta)
  eta.vec <- jj$arg1
  attr <- jj$attr
  
  eta.out <- rep(eta.vec,(kmax+1))
  jj <- .C("coulomb_CL_array",
           as.double(L_min),
           as.integer(kmax),
           as.double(eta.vec),
           as.integer(length(eta.vec)),
           val=as.double(eta.out),
           status=as.integer(0*eta.vec),
           PACKAGE="gsl"
           )

  val <- jj$val
  dim(val) <- c(kmax+1, length(eta.vec))
  rownames(val) <- L_min:(L_min+kmax)
  status <- jj$status
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }
  
  if(give){
    return(list(val=val, status=status))
  } else {
    return(val)
  }  
}  
