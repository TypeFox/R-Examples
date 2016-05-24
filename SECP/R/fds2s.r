# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# fds2s() and fds3s() functions use a linear regression model for
# statistical estimation of the mass fractal dimension of sampling
# clusters on 2D & 3D square lattice. 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# rfq - matrix of relative sampling frequencies 
#       for sites of the 2D & 3D percolation lattice;
# bnd - list of boundary coordinates and radii 
#       for the isotropic 2D & 3D set cover.
# Variables:
# w, r - relative frequency sums and lengths for 
#        the iso- & anisotropic 2D & 3D set cover.
# Value:
# - linear regression model for statistical estimation of the mass
#   fractal dimension of sampling clusters on 2D & 3D square lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
require(SPSL)
fds2s <- function(rfq=fssi20(x=95), 
                  bnd=isc2s(k=12, x=dim(rfq))) {
  w <- rep(0, times=ncol(bnd))
  for (i in seq(ncol(bnd)))
    w[i] <- sum(rfq[bnd["x1",i]:bnd["x2",i],
                    bnd["y1",i]:bnd["y2",i]])
  r <- log(bnd["r",])
  w <- log(w)
  return(lm(w ~ r))
}
fds3s <- function(rfq=fssi30(x=95),
                  bnd=isc3s(k=12, x=dim(rfq))) {
  w <- rep(0, times=ncol(bnd))
  for (i in seq(ncol(bnd)))
    w[i] <- sum(rfq[bnd["x1",i]:bnd["x2",i],
                    bnd["y1",i]:bnd["y2",i],
                    bnd["z1",i]:bnd["z2",i]])
  r <- log(bnd["r",])
  w <- log(w)
  return(lm(w ~ r))
}