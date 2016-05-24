require(fields)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function:
# voss1g() & voss2g() - a fractal Brownian function on uniform
# 1D & 2D grid with a generic version of the Voss algorithm.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# p - a matrix of parameters: 
# nrow(p) - a number of iterations;
# p[,"n"] - a number of partition points in the iteration process;
# p[,"s"] - a standard deviation of normal pseudorandom additions;
# center - logical; if center=TRUE, y-coordinates of
# prefractal points will be centered.
# Value:
# voss - a list of Cartesian coordinates of prefractal points.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
voss1g <- function(p=cbind(n=0.5^-seq(0,7)+1, 
                                s=dchisq(seq(0,7), df=2)), 
                        center=TRUE) {
  voss <- list(x=c(0,1), y=c(0,0))
  for (i in seq(2,nrow(p))) { 
    voss <- approx(voss, n=p[i,"n"])
    voss$y <- rnorm(n=p[i,"n"], mean=voss$y, sd=p[i,"s"])
  }
  if (center) voss$y <- voss$y - mean(voss$y)
  return(voss)
}
voss2g <- function(p=cbind(n=0.5^-seq(0,7)+1, 
                                s=dchisq(seq(0,7), df=2)), 
                        center=TRUE) {
  voss <- list(x=c(0,1), y=c(0,1), z=array(0, dim=c(2,2)))
  for (i in seq(2,nrow(p))) { 
    gxy <- seq(0, 1, length=p[i,"n"])
    voss <- interp.surface.grid(voss, list(x=gxy, y=gxy))
    voss$z <- array(rnorm(n=p[i,"n"]^2, mean=voss$z, sd=p[i,"s"]), 
                    dim=rep(p[i,"n"], times=2)) 
  }
  if (center) voss$z <- voss$z - mean(voss$z)
  return(voss) 
}