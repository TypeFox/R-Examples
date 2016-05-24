CompareL1_L2_TV1 <- function(N,M,per) {
  # Part of R1Magic by mehmet.suzen@physics.org
  phi <- matrix( rnorm(N * M ),  M,N )
  xorg <-  sparseSignal(N, per*N)
  y <- phi %*% xorg
  T <- diag(N)
  p <- matrix(0,N,1)
  ll<-solveL1(phi,y,T,p)
  x1<-ll$estimate
  ll<-solveL2(phi,y,T,p)
  x2<-ll$estimate
  ll<-solve1TV(phi,y,T,p)
  xtv<-ll$estimate
 return( data.frame(xorg,x1,x2,xtv) )
}

