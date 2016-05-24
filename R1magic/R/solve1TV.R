solve1TV <- function(phi,y,T,x0,lambda=0.1) {
  # Part of R1Magic by mehmet.suzen@physics.org
  return( nlm(objective1TV, p=x0, T=T, phi=phi, y=y, lambda=lambda ) )
}

