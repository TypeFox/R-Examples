ph2single <- function(pu,pa,ep1,ep2,nsoln=5) {
  n0 <- 1
  incr.n0 <- TRUE
  qep1 <- qbinom(1-ep1, n0, pu)
  err2 <- pbinom(qep1, n0, pa)
  if (err2 <= ep2) incr.n0 <- FALSE
  while (incr.n0) {
    n0 <- n0 + 1
    qep1 <- qbinom(1-ep1, n0, pu)
    err2 <- pbinom(qep1, n0, pa)
    if (err2 <= ep2) incr.n0 <- FALSE
  }
  n <- n0
  soln <- matrix(0,ncol=4,nrow=nsoln)
  isoln <- 0
  while(isoln < nsoln) {
    for(r in 0:(n-1)) {
      err1 <- 1-pbinom(r,n,pu) 
      err2 <- pbinom(r,n,pa)
      if((err1 <= ep1) & (err2 <= ep2)) {
        isoln <- isoln + 1
        if(isoln <= nsoln) soln[isoln,] <- c(n,r,err1,err2)
      }
    }
    n <- n + 1
  }
  soln <- as.data.frame(soln)
  names(soln) <- c("n","r","Type I error","Type II error")
  soln
}
