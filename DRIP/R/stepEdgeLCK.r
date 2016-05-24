# This is R source code for function 'stepEdgeLCK', in the R
# package "image".
# Date: April 25, 2013
# Creator: Yicheng Kang

stepEdgeLCK = function(image, bandwidth, thresh, plot =
  FALSE){
  if (!is.matrix(image))
    stop("image data must be a matrix")
  n1 = as.integer(dim(image)[1])
  n2 = as.integer(dim(image)[2])
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
  if (as.integer(bandwidth) < 1)
    stop("bandwidth must be a positive integer")
  if (n1 + 2 * bandwidth + 2 > 600)
    stop("bandwidth is too large or the resolution
of the image is too high.")
  if (!is.numeric(thresh))
    stop("threshold  must be numeric")
  n1 = dim(image)[1]
  z = matrix(as.double(image), ncol = n1)
  k = as.integer(bandwidth)
  u = as.double(thresh)
  out = .Fortran('lck_diff', n = as.integer(n1 - 1), obsImg = z,
    bandwidth = as.integer(k), diff = z)
  edge = out$diff
  edge[out$diff < u] = 0
  edge[out$diff >= u] = 1
  if (plot == FALSE) { return(edge) }
  if (plot == TRUE) {
    x = seq(0, 1, length = n1); y = x
    image(x, y, 1 - edge, col = gray(0:1))
    return(edge)
  }
}
