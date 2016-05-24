require(fields)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function:
# voss1d() & voss2d() - a fractal Brownian function on uniform 
# 1D & 2D grid with a classic version of the Voss algorithm.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# g - a number of iteration; H - a Hurst parameter; 
# r - a partition coefficient for the iteration segment; 
# center - logical; if center=TRUE, y-coordinates of
# prefractal points will be centered.
# Variables:
# n - a number of partition points in the iteration process;
# s - a standard deviation of normal pseudorandom additions.
# Value:
# voss - a list of Cartesian coordinates of prefractal points.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
voss1d <- function(g=7, H=0.5, r=0.5, center=TRUE) {
  n <- r^-seq(0,g)+1
  s <- H*log(1/r)*r^(H*seq(0,g))
  voss <- list(x=c(0,1), y=c(0,0))
  for (i in seq(2,g+1)) {
    voss <- approx(voss, n=n[i])
    voss$y <- rnorm(n=n[i], mean=voss$y, sd=s[i])
  }
  if (center) voss$y <- voss$y - mean(voss$y)
  return(voss)
}
voss2d <- function(g=7, H=0.5, r=0.5, center=TRUE) {
  n <- r^-seq(0,g)+1
  s <- H*log(1/r)*r^(H*seq(0,g))
  voss <- list(x=c(0,1), y=c(0,1), z=array(0, dim=c(2,2)))
  for (i in seq(2,g+1)) {
    gxy <- seq(0, 1, length=n[i])
    voss <- interp.surface.grid(voss, list(x=gxy, y=gxy))
    voss$z <- array(rnorm(n=n[i]^2, mean=voss$z, sd=s[i]), 
                    dim=rep(n[i], times=2)) 
  }
  if (center) voss$z <- voss$z - mean(voss$z)
  return(voss) 
}