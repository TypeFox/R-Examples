# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# License: GPL-3
# Package: Site Percolation on Square Lattice (SPSL)
# Author: Pavel V. Moskalev <moskalefff@gmail.com>
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# ssi2d() and ssi3d() functions provide a labeling of 
# isotropic 2D & 3D clusters with Moore d-neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# x - linear dimension of the percolation lattice; 
# p0 - relative fraction of accessible sites
#      (occupation probability) for percolation lattice;
# p1, p2 - p1 value, weighted by 2D & 3D d-neighborhood;
# set - vector of linear indexes of initial sites subset;
# all - trigger "Mark all initial sites or accessible only?"
# Variables:
# e0, e1, e2 - linear indexes of sites combinations from 
#              2D & 3D Moore neighborhood;
# b - length of initial sites subset.
# Value:
# acc - labeled accessibility matrix for the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
ssi2d <- function(x=33, 
                  p0=0.5, p1=p0/2, 
                  set=(x^2+1)/2, all=TRUE) {
  e0 <- c(-1, 1,-x, x)
  e1 <- colSums(matrix(e0[c(
    1,3, 2,3, 1,4, 2,4)], nrow=2))
  e <- as.integer(c(e0,e1))
  p0 <- rep(p0, length(e0))
  p1 <- rep(p1, length(e1))
  p <- as.double(c(p0,p1))
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
ssi3d <- function(x=33, 
                  p0=0.2, p1=p0/2, p2=p0/3,
                  set=(x^3+1)/2, all=TRUE) {
  e0 <- c(-1, 1,-x, x,-x^2, x^2)
  e1 <- colSums(matrix(e0[c(
    1,3, 2,3, 1,4, 2,4,
    1,5, 2,5, 1,6, 2,6,
    3,5, 4,5, 3,6, 4,6)], nrow=2))
  e2 <- colSums(matrix(e0[c(
    1,3,5, 2,3,5, 1,4,5, 2,4,5,
    1,3,6, 2,3,6, 1,4,6, 2,4,6)], nrow=3))
  e <- as.integer(c(e0,e1,e2))
  p0 <- rep(p0, length(e0))
  p1 <- rep(p1, length(e1))
  p2 <- rep(p2, length(e2))
  p <- as.double(c(p0,p1,p2))  
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