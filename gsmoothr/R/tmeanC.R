tmeanC <- function(sp, x, spout=NULL, nProbes=10, probeWindow=600, trim=0.1) {
  # void tmean(double *x, double *xs, long *sp, long *n_, double *tr, long *_np, long *_pw)

  stopifnot( length(x)==length(sp) )
  
  k <- !is.na(x) & !is.nan(x)
  ni <- sum(k)

  if( is.null(spout) ) {
    outx <- rep(NA,ni)
    # screen out NA and NaNs before sending to C
    out<-.C("tmean", x=as.double(x[k]), xs=double(ni), sp=as.integer(sp[k]), n=as.integer(ni), 
                     tr=as.double(trim), np=as.integer(nProbes), pw=as.integer(probeWindow), PACKAGE="gsmoothr")
    outx[k] <- out$xs
  } else {
    no <- length(spout)
    outx <- rep(NA,no)
    out<-.C("tmeanPos", xi=as.double(x[k]), xo=double(no), spi=as.integer(sp[k]), spo=as.integer(spout), 
                        ni=as.integer(ni), no=as.integer(no), tr=as.double(trim), np=as.integer(nProbes), 
                        pw=as.integer(probeWindow), PACKAGE="gsmoothr")
    outx <- out$xo
  }

  outx
}
