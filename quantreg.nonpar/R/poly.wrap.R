poly.wrap <-
function(x, degree=1, coefs=NULL, nderivs=1, raw=FALSE)
{
  if (nderivs==0){
    out <- poly(x=x, degree=degree, coefs=coefs, raw=raw)
  } else if (nderivs==1) {
    out <- dpoly(x=x, degree=degree, coefs=coefs, raw=raw)
  } else if (nderivs==2) {
    out <- ddpoly(x=x, degree=degree, coefs=coefs, raw=raw)
  } else {
    stop("This method is only available for nderivs=0, nderivs=1, and nderivs=2")
  }

  return(out)
}
