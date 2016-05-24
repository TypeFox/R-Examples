# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# fdc2s() and fdc3s() functions use a linear regression model for 
# statistical estimation of the mass fractal dimension of a site 
# cluster on 2D & 3D square lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# acc - labeled matrix for sites of 
#       the 2D & 3D percolation lattice.
# bnd - list of boundary coordinates and radii 
#       for the isotropic 2D & 3D set cover.
# Variables:
# n, r - absolute frequency sums and lengths for 
#        the iso- & anisotropic 2D & 3D set cover.
# Value:
# - linear regression model for statistical estimation of the mass
#   fractal dimension of a site cluster on 2D & 3D square lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
require(SPSL)
fdc2s <- function(acc=ssi20(x=95), 
                  bnd=isc2s(k=12, x=dim(acc))) {
  n <- rep(0, times=ncol(bnd))
  for (i in seq(ncol(bnd)))
    n[i] <- sum(acc[bnd["x1",i]:bnd["x2",i],
                    bnd["y1",i]:bnd["y2",i]] > 1)
  r <- log(bnd["r",])
  n <- log(n)
  return(lm(n ~ r))
}
fdc3s <- function(acc=ssi30(x=95),
                  bnd=isc3s(k=12, x=dim(acc))) {
  n <- rep(0, times=ncol(bnd))
  for (i in seq(ncol(bnd)))
    n[i] <- sum(acc[bnd["x1",i]:bnd["x2",i],
                    bnd["y1",i]:bnd["y2",i],
                    bnd["z1",i]:bnd["z2",i]] > 1)
  r <- log(bnd["r",])
  n <- log(n)
  return(lm(n ~ r))
}