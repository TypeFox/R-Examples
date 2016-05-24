xewma.q.crit.prerun <- function(l, L0, mu, p, zr=0, hs=0, sided="two", limits="fix", size=100, df=NULL, estimated="mu", qm.mu=30, qm.sigma=30,
                                truncate=1e-10, bound=1e-10, c.error=1e-10, p.error=1e-9, OUTPUT=FALSE) {                                
                                
  if ( OUTPUT ) cat("\nc\t\tp\n")
  
  c2 <- xewma.q.crit(l, L0, mu, p, zr=zr, hs=hs, sided=sided, limits=limits, OUTPUT=FALSE)  
  p2 <- 1 - xewma.sf.prerun(l, c2, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1,
                            size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate, bound=bound)[L0]                            
  if ( OUTPUT ) cat(paste(c2,"\t",p2,"\n"))
  
  if ( p2 > p ) {
    while ( p2 > p ) {
      p1 <- p2
      c2 <- c2 + .5
      p2 <- 1 - xewma.sf.prerun(l, c2, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1,
                                size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate, bound=bound)[L0]
      if ( OUTPUT ) cat(paste(c2,"\t",p2,"\n"))
    }
    c1 <- c2 - .5
  } else {
    while ( p2 <= p ) {
      p1 <- p2
      c2 <- c2 - .5
      p2 <- 1 - xewma.sf.prerun(l, c2, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1,
                                size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate, bound=bound)[L0]
      if ( OUTPUT ) cat(paste(c2,"\t",p2,"\n"))
    }
    c1 <- c2 + .5
  }
  
  if ( size < 41 ) {
    if ( qm.mu < 70 ) qm.mu <- 70 
    if ( qm.mu < 70 ) qm.mu <- 70
    if ( size < 21 ) {
      if ( qm.mu < 90 ) qm.mu <- 90 
      if ( qm.mu < 90 ) qm.mu <- 90
    }
    if ( p2 > p ) {
      while ( p2 > p ) {
        p1 <- p2
        c2 <- c2 + .1
        p2 <- 1 - xewma.sf.prerun(l, c2, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1,
                                  size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate, bound=bound)[L0]
        if ( OUTPUT ) cat(paste(c2,"\t",p2,"\n"))
      }
      c1 <- c2 - .1
    } else {
      while ( p2 <= p ) {
        p1 <- p2
        c2 <- c2 - .1
        p2 <- 1 - xewma.sf.prerun(l, c2, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1,
                                  size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate, bound=bound)[L0]
        if ( OUTPUT ) cat(paste(c2,"\t",p2,"\n"))
      }
      c1 <- c2 + .1
    }
  }
  
  p.error_ <- 1; c.error_ <- 1
  while ( p.error_ > p.error & c.error_ > c.error ) {
    c3 <- c1 + (p - p1)/(p2 - p1)*(c2 - c1)
    p3 <- 1 - xewma.sf.prerun(l, c3, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1,
                              size=size, df=df, estimated=estimated, qm.mu=qm.mu, qm.sigma=qm.sigma, truncate=truncate, bound=bound)[L0]
    if ( OUTPUT ) cat(paste(c3,"\t",p3,"\n"))
    c1 <- c2; c2 <- c3
    p1 <- p2; p2 <- p3
    p.error_ <- abs(p2 - p); c.error_ <- abs(c2 - c1)
  }
  
  names(c3) <- "c"
  c3
}
