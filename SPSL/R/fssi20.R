# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# License: GPL-3
# Package: Site Percolation on Square Lattice (SPSL)
# Author: Pavel V. Moskalev <moskalefff@gmail.com>
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# fssi20() and fssi30() functions calculates the relative 
# frequency distribution of isotropic 2D & 3D clusters 
# with von Neumann neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# n - sample size;
# x - linear dimension of the percolation lattice; 
# p - relative fraction of accessible sites 
#     (occupation probability) for percolation lattice;
# set - vector of linear indexes of initial sites subset;
# all - trigger "Mark all initial sites or accessible only?"
# Value:
# rfq - matrix of relative sampling frequencies for sites of
#       the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fssi20 <- function(n=1000, 
                   x=33, p=0.592746, 
                   set=(x^2+1)/2, all=TRUE) {
  rfq <- array(0, dim=rep(x, times=2))
  for (i in seq(n))
    rfq <- rfq + (ssi20(x, p, set, all) > 1)
  return(rfq/n)
}
fssi30 <- function(n=1000,
                   x=33, p=0.311608, 
                   set=(x^3+1)/2, all=TRUE) {
  rfq <- array(0, dim=rep(x, times=3))
  for (i in seq(n))
    rfq <- rfq + (ssi30(x, p, set, all) > 1)
  return(rfq/n)
}