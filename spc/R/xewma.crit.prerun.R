xewma.crit.prerun <- function(l, L0, mu, zr=0, hs=0, sided="two", limits="fix", size=100, df=NULL, estimated="mu", qm.mu=30, qm.sigma=30,
                              truncate=1e-10, c.error=1e-12, L.error=1e-9, OUTPUT=FALSE) {
                              
  if ( OUTPUT ) cat("\nc\t\tL\n")
  
  c2 <- xewma.crit(l, L0, mu0=mu, zr=zr, hs=hs, sided=sided, limits=limits)
  L2 <- xewma.arl.prerun(l, c2, mu, zr=zr, hs=hs, sided=sided, limits=limits, q=1, size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate)
  if ( OUTPUT ) cat(paste(c2,"\t",L2,"\n"))
  
  if ( L2 < L0 ) {
    while ( L2 < L0 ) {
      L1 <- L2
      c2 <- c2 + .5
      L2 <- xewma.arl.prerun(l, c2, mu, zr=zr, hs=hs, sided=sided, limits=limits, q=1, size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate)
      if ( OUTPUT ) cat(paste(c2,"\t",L2,"\n"))
    }
    c1 <- c2 - .5
  } else {
    while ( L2 >= L0 ) {
      L1 <- L2
      c2 <- c2 - .5
      L2 <- xewma.arl.prerun(l, c2, mu, zr=zr, hs=hs, sided=sided, limits=limits, q=1, size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate)
      if ( OUTPUT ) cat(paste(c2,"\t",L2,"\n"))
    }
    c1 <- c2 + .5
  }
  
  if ( size < 51 ) {
    if ( qm.mu < 70 ) qm.mu <- 70 
    if ( qm.mu < 70 ) qm.mu <- 70
    if ( size < 31 ) {
      if ( qm.mu < 90 ) qm.mu <- 90 
      if ( qm.mu < 90 ) qm.mu <- 90
    }
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        c2 <- c2 + .1
        L2 <- xewma.arl.prerun(l, c2, mu, zr=zr, hs=hs, sided=sided, limits=limits, q=1, size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate)
        if ( OUTPUT ) cat(paste(c2,"\t",L2,"\n"))
      }
      c1 <- c2 - .1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        c2 <- c2 - .1
        L2 <- xewma.arl.prerun(l, c2, mu, zr=zr, hs=hs, sided=sided, limits=limits, q=1, size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate)
        if ( OUTPUT ) cat(paste(c2,"\t",L2,"\n"))
      }
      c1 <- c2 + .1
    }
  }

  L.error_ <- 1; c.error_ <- 1
  while ( L.error_ > L.error & c.error_ > c.error ) {
    c3 <- c1 + (L0 - L1)/(L2 - L1)*(c2 - c1)
    L3 <- xewma.arl.prerun(l, c3, mu, zr=zr, hs=hs, sided=sided, limits=limits, q=1, size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate)
    if ( OUTPUT ) cat(paste(c3,"\t",L3,"\n"))
    c1 <- c2; c2 <- c3
    L1 <- L2; L2 <- L3
    L.error_ <- abs(L2 - L0); c.error_ <- abs(c2 - c1)
  }
  
  names(c3) <- "c"
  c3
}
