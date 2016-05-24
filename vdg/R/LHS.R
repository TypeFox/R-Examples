#' Latin Hypercube Sampling
#' 
#' Different versions of latin hypercube sampling (LHS): ordinary LHS,  midpoint LHS, symmetric LHS or randomized symmetric LHS. LHS is a method 
#' for constructing space-filling designs. They can be more efficient for hypercuboidal design regions than other sampling methods.
#' 
#' @aliases LHS MLHS SLHS RSLHS
#' @param n number of design points to generate
#' @param m number of design factors
#' @param lim limits of the coordinates in all dimensions
#' @return Matrix with samples as rows. 
#' @author Pieter C. Schoonees
#' @examples
#' 
#' set.seed(1234)
#' pts <- seq(-1, 1, length = 11)
#' 
#' # Ordinary LHS
#' samp <- LHS(n = 10, m = 2)
#' plot(samp, main = "LHS")
#' abline(h = pts, v = pts, col = "lightgrey")
#' 
#' # Midpoint LHS
#' samp <- MLHS(n = 10, m = 2)
#' plot(samp, main = "MLHS")
#' abline(h = pts, v = pts, col = "lightgrey")
#' 
#' # Symmetric LHS
#' samp <- SLHS(n = 10, m = 2)
#' plot(samp, main = "SLHS")
#' abline(h = pts, v = pts, col = "lightgrey")
#' 
#' # Randomized Symmetric LHS
#' samp <- RSLHS(n = 10, m = 2)
#' plot(samp, main = "RSLHS")
#' abline(h = pts, v = pts, col = "lightgrey")
#' @keywords design
#' @export
LHS <-
function (n, m = 3, lim = c(-1, 1)) 
{
  pts <- seq(from = lim[1], to = lim[2], length = n + 1)
  pts <- pts[-1]
  samp <- matrix(0, nrow = n, ncol = m)
  for(i in 1:m)
  	samp[,i] <- pts[sample(1:n, n)]
  umat <- matrix(runif(m*n), nrow = n, ncol = m)
  samp <- samp - (lim[2] - lim[1])*umat/n
  samp
}
#' @rdname LHS
#' @export
MLHS <-
  function (n, m = 3, lim = c(-1, 1)) 
  {
    pts <- seq(from = lim[1], to = lim[2], length = n + 1)
    pts <- pts[-1]
    samp <- matrix(0, nrow = n, ncol = m)
    for(i in 1:m)
      samp[,i] <- pts[sample(1:n, n)]
    samp <- samp - (lim[2] - lim[1])*0.5/n
    samp
  }
#' @rdname LHS
#' @export
SLHS <-
  function (n, m = 3, lim = c(-1, 1)) 
  {
    k <- n/2
    if(n %% 2 != 0) stop("Not a even number of points - a symmetric LHD cannot be constructed.")
    pts <- seq(from = lim[1], to = lim[2], length = n + 1)
    pts <- pts[-1]
    samp <- matrix(0, nrow = k, ncol = m)
    for(i in 1:m)
      samp[,i] <- sample(1:n, k)	
    samp <- rbind(samp, n + 1 - samp)
    for(i in 1:m)
      samp[,i] <- pts[samp[,i]]
    samp <- samp - (lim[2] - lim[1])*0.5/n
    samp
  }
#' @rdname LHS
#' @export
RSLHS <-
  function (n, m = 3, lim = c(-1, 1)) 
  {
    k <- n/2
    if(n %% 2 != 0) stop("Not a even number of points - a symmetric LHD cannot be constructed.")
    pts <- seq(from = lim[1], to = lim[2], length = n + 1)
    pts <- pts[-1]
    samp <- matrix(0, nrow = k, ncol = m)
    for(i in 1:m)
      samp[,i] <- sample(1:n, k)	
    samp <- rbind(samp, n + 1 - samp)
    # samp2 <- samp/(n/2)-1
    for(i in 1:m)
      samp[,i] <- pts[samp[,i]]
    umat <- matrix(runif(n*m), nrow = n, ncol = m)
    samp <- samp - (lim[2] - lim[1])*umat/n
    samp
  }

