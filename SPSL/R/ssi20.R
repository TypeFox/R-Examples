# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# License: GPL-3
# Package: Site Percolation on Square Lattice (SPSL)
# Author: Pavel V. Moskalev <moskalefff@gmail.com>
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# ssi20() and ssi30() functions provide a labeling of 
# isotropic 2D & 3D clusters with von Neumann neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# x - linear dimension of the percolation lattice; 
# p - relative fraction of accessible sites 
#     (occupation probability) for percolation lattice;
# set - vector of linear indexes of initial sites subset;
# all - trigger "Mark all initial sites or accessible only?"
# Variables:
# e - linear indexes of sites from 2D & 3D von Neumann neighborhood;
# b - length of initial sites subset.
# Value:
# acc - labeled accessibility matrix for the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
ssi20 <- function(x=33, 
                  p=0.592746, 
                  set=(x^2+1)/2, all=TRUE) {
  e <- as.integer(c(-1, 1,-x, x))
  p <- as.double(rep(p, length(e)))
  acc <- array(runif(x^2), rep(x,2))
  if (!all) set <- set[acc[set] < mean(p)]
  b <- as.integer(length(set))
  cls <- rep(0L, max(p)*x^2 + b*all)
  acc[set] <- 2 
  acc[c(1,x),] <- acc[,c(1,x)] <- 1
  cls[seq_along(set)] <- as.integer(set - 1)
  .Call("ssTNd", p, acc, b, e, cls) 
  return(acc)
}
ssi30 <- function(x=33, 
                  p=0.311608, 
                  set=(x^3+1)/2, all=TRUE) {
  e <- as.integer(c(-1, 1,-x, x,-x^2, x^2))
  p <- as.double(rep(p, length(e)))
  acc <- array(runif(x^3), rep(x,3))
  if (!all) set <- set[acc[set] < mean(p)]
  b <- as.integer(length(set))
  cls <- rep(0L, max(p)*x^3 + b*all)
  acc[set] <- 2 
  acc[c(1,x),,] <- acc[,c(1,x),] <- acc[,,c(1,x)] <- 1
  cls[seq_along(set)] <- as.integer(set - 1)
  .Call("ssTNd", p, acc, b, e, cls)   
  return(acc)
}