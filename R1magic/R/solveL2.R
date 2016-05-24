solveL2 <- function(phi,y,T,x0,lambda=0.1) {
  # Part of R1Magic by mehmet.suzen@physics.org
  return( nlm(objectiveL2, p=x0, T=T, phi=phi, y=y, lambda=lambda ) )
}

