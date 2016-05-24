# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# License: GPL-3
# Package: Site Percolation on Square Lattice (SPSL)
# Author: Pavel V. Moskalev <moskalefff@gmail.com>
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# fssi2d() and fssi3d() functions calculates the relative 
# frequency distribution of isotropic 2D & 3D clusters 
# with Moore d-neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# n - sample size;
# x - linear dimension of the percolation lattice; 
# p0 - relative fraction of accessible sites
#      (occupation probability) for percolation lattice;
# p1, p2 - p1 value, weighted by 2D & 3D d-neighborhood;
# set - vector of linear indexes of initial sites subset;
# all - trigger "Mark all initial sites or accessible only?"
# Value:
# rfq - matrix of relative frequencies for sites of the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fssi2d <- function(n=1000, 
                   x=33, p0=0.5, p1=p0/2,
                   set=(x^2+1)/2, all=TRUE) {
  rfq <- array(0, dim=rep(x, times=2))
  for (i in seq(n))
    rfq <- rfq + (ssi2d(x, p0, p1, set, all) > 1)
  return(rfq/n)
}
fssi3d <- function(n=1000,
                   x=33, p0=0.2, p1=p0/2, p2=p0/3, 
                   set=(x^3+1)/2, all=TRUE) {
  rfq <- array(0, dim=rep(x, times=3))
  for (i in seq(n))
    rfq <- rfq + (ssi3d(x, p0, p1, p2, set, all) > 1)
  return(rfq/n)
}