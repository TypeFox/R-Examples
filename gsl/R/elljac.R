"elljac" <- function(u, m, give=FALSE, strict=TRUE){
  jj <- process.args(u, m)
  u.vec <- jj$arg1
  m.vec <- jj$arg2
  attr <- jj$attr

  jj <- .C("elljac_e",
           as.double(u.vec),
           as.double(m.vec),
           as.integer(length(u.vec)),
           sn=as.double(u.vec),
           cn=as.double(u.vec),
           dn=as.double(u.vec),
           status=as.integer(0*u.vec),
           PACKAGE="gsl"
           )

  sn <- jj$sn
  cn <- jj$cn
  dn <- jj$dn
  attributes(sn) <- attr
  attributes(cn) <- attr
  attributes(dn) <- attr

  status <- jj$status
  attributes(status) <- attr
  
  if(strict){
    sn <- strictify(sn,status)
    cn <- strictify(cn,status)
    dn <- strictify(dn,status)
  }
  
  if(give){
      return(list(list(sn=sn,cn=cn,dn=dn),status=status))
  } else {
    return(list(sn=sn,cn=cn,dn=dn))
  }
}

"sn_cn_dn" <- function(z,m,thing){  # complex case
  jj.r <- elljac(Re(z),m)
  s <- jj.r$sn
  c <- jj.r$cn
  d <- jj.r$dn
  if(is.complex(z)){
    jj.i <- elljac(Im(z),1-m)
    s1 <- jj.i$sn
    c1 <- jj.i$cn
    d1 <- jj.i$dn
    out <- switch(thing,
                  sn = (s*d1   +1i*c*d*s1*c1)/(c1^2+m*s^2*s1^2),
                  cn = (c*c1   -1i*s*d*s1*d1)/(c1^2+m*s^2*s1^2),
                  dn = (d*c1*d1-1i*m*s*c*s1 )/(c1^2+m*s^2*s1^2),
                  stop('argument "thing" should be one of sn, cn, dn')
                  )
  } else {
    out <- switch(thing,
                  sn = s,
                  cn = c,
                  dn = d,
                  stop('argument "thing" should be one of sn, cn, dn')
                  )
  }
  return(out)
}

gsl_sn <- function(z,m){sn_cn_dn(z,m,thing="sn")}
gsl_cn <- function(z,m){sn_cn_dn(z,m,thing="cn")}
gsl_dn <- function(z,m){sn_cn_dn(z,m,thing="dn")}

gsl_ns <- function(z,m){1/gsl_sn(z,m)}
gsl_nc <- function(z,m){1/gsl_cn(z,m)}
gsl_nd <- function(z,m){1/gsl_dn(z,m)}

gsl_sc <- function(z,m){gsl_sn(z,m)/gsl_cn(z,m)}
gsl_sd <- function(z,m){gsl_sn(z,m)/gsl_dn(z,m)}

gsl_cs <- function(z,m){gsl_cn(z,m)/gsl_sn(z,m)}
gsl_cd <- function(z,m){gsl_cn(z,m)/gsl_dn(z,m)}

gsl_ds <- function(z,m){gsl_dn(z,m)/gsl_sn(z,m)}
gsl_dc <- function(z,m){gsl_dn(z,m)/gsl_cn(z,m)}
