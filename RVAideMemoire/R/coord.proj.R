coord.proj <- function (coord,slp) {
  coeff <- matrix(c(-slp,1/slp,1,1),ncol=2,dimnames=list(1:2,c("x","y")))
  d <- integer(nrow(coord))
  for (i in 1:nrow(coord)) {
    int.perp <- coord[i,2]+1/slp*coord[i,1]
    ints <- matrix(c(0,int.perp),ncol=1)
    val <- solve(coeff,ints)
    xproj <- val[1,1]
    yproj <- val[2,1]
    d[i] <- sqrt(xproj^2+yproj^2)*sign(xproj)
  }
  return(d)
}
