"coef.bigqr" <-
function(bigQR,nvar=NULL,...){
  p <- length(bigQR$D)
  if (is.null(nvar))
    nvar <- p
  if (nvar <1 | nvar >p) stop("Invalid value of `nvar'")

  if (!bigQR$checked)
    bigQR<-singcheck.bigqr(bigQR)
  
  tmp <- .Fortran("regcf", as.integer(p),
                  as.integer(p*p/2),
                  bigQR$D,
                  bigQR$rbar,
                  bigQR$thetab,
                  bigQR$tol,
                  beta=numeric(p),
                  nreq=as.integer(nvar),
                  ier=integer(1), DUP=FALSE)

  if (tmp$ier!=0) stop("Error in REGCF: can't happen")

  tmp$beta
}
