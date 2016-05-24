xewma.q.crit <- function(l, L0, mu, alpha, zr=0, hs=0, sided="two", limits="fix", r=40, c.error=1e-12, a.error=1e-9, OUTPUT=FALSE) {

  c2 <- 0
  p2 <- 1
  if ( OUTPUT ) cat("\nc\t\tp\n")
  while ( p2 > alpha ) {
    p1 <- p2
    c2 <- c2 + .5
    p2 <- 1 - xewma.sf(l, c2, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1, r=r)[L0]
    if ( OUTPUT ) cat(paste(c2,"\t",p2,"\n"))
  }
  while ( p2 <= alpha & c2 > .02 ) {
    p1 <- p2
    c2 <- c2 - .02
    p2 <- 1 - xewma.sf(l, c2, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1, r=r)[L0]
    if ( OUTPUT ) cat(paste(c2,"\t",p2,"\n"))
  }
  c1 <- c2 + .02
  
  a.error_ <- 1; c.error_ <- 1
  while ( a.error_ > a.error & c.error_ > c.error ) {
    c3 <- c1 + (alpha - p1)/(p2 - p1)*(c2 - c1)
    p3 <- 1 - xewma.sf(l, c3, mu, L0, zr=zr, hs=hs, sided=sided, limits=limits, q=1, r=r)[L0]
    if ( OUTPUT ) cat(paste(c3,"\t",p3,"\n"))
    c1 <- c2; c2 <- c3
    p1 <- p2; p2 <- p3
    a.error_ <- abs(p2 - alpha); c.error_ <- abs(c2 - c1)
  }
  
  names(c3) <- "c"
  c3
}
