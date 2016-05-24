
xDshewhartrunsrules.arl <- function(delta, c=1, m=NULL, type="12") {
  eps <- 1e-6
  if ( is.null(m) ) {
    m <- 4
    arl1 <- xDshewhartrunsrulesFixedm.arl(delta, c=c, m=m, type=type)
    arl2 <- arl1 + 2*eps
    while ( abs(arl2-arl1)>eps & m<1e4 ) {
      m <- round(1.5 * m)
      arl1 <- xDshewhartrunsrulesFixedm.arl(delta, c=c, m=m, type=type)
      arl2 <- xDshewhartrunsrulesFixedm.arl(delta, c=c, m=m+1, type=type)
    }
    arl <- arl1
  } else {
    arl <- xDshewhartrunsrulesFixedm.arl(delta, c=c, m=m, type=type)
  }
  arl
}