
xshewhartrunsrules.crit <- function(L0, mu=0, type="12") {	
  if (type=="14" & L0>255) {
    stop("L0 too large for type=\"14\"")
  } else {
    c1 <- 1
    c2 <- 1.5
    arl1 <- xshewhartrunsrules.arl(mu,c=c1,type=type)
    arl2 <- xshewhartrunsrules.arl(mu,c=c2,type=type)
    a.error <- 1; c.error <- 1
    while ( a.error>1e-6 && c.error>1e-8 ) {
      c3 <- c1 + (L0-arl1)/(arl2-arl1)*(c2-c1)
      arl3 <- xshewhartrunsrules.arl(mu,c=c3,type=type)
      c1 <- c2; c2 <- c3
      arl1 <- arl2; arl2 <- arl3
      a.error <- abs(arl2-L0); c.error <- abs(c2-c1)
    }
  }
  c3
}