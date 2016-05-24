# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# License: GPL-3
# Package: Site Percolation on Square Lattice (SPSL)
# Author: Pavel V. Moskalev <moskalefff@gmail.com>
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# ssa2d() and ssa3d() functions provide a labeling of 
# anisotropic 2D & 3D clusters with d-neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# x - linear dimension of the percolation lattice; 
# p0 - vector of p-values, distributed by lattice 
#      directions: (-x, +x, -y, +y, -z, +z);
# p1, p2 - double and triple combinations of p0-components, 
#          weighted by 2D & 3D d-neighborhood;
# set - vector of linear indexes of initial sites subset;
# all - trigger "Mark all initial sites or accessible only?"
# Variables:
# e0, e1, e2 - linear indexes of sites combinations from 
#              2D & 3D Moore neighborhood;
# b - length of initial sites subset.
# Value:
# acc - labeled accessibility matrix for the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
ssa2d <- function(x=33, 
                  p0=runif(4, max=0.8), 
                  p1=colMeans(matrix(p0[c(
                    1,3, 2,3, 1,4, 2,4)], nrow=2))/2,
                  set=(x^2+1)/2, all=TRUE) {
  e0 <- c(-1, 1,-x, x)
  e1 <- colSums(matrix(e0[c(
    1,3, 2,3, 1,4, 2,4)], nrow=2))
  e <- as.integer(c(e0,e1))
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
ssa3d <- function(x=33, 
                  p0=runif(6, max=0.4),
                  p1=colMeans(matrix(p0[c(
                    1,3, 2,3, 1,4, 2,4,
                    1,5, 2,5, 1,6, 2,6,
                    3,5, 4,5, 3,6, 4,6)], nrow=2))/2,
                  p2=colMeans(matrix(p0[c(
                    1,3,5, 2,3,5, 1,4,5, 2,4,5,
                    1,3,6, 2,3,6, 1,4,6, 2,4,6)], nrow=3))/3,                  
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