# This is R source code for function 'diffLL2K' in the R
# package 'image'.
# Creator: Yicheng Kang
# Date: April 29, 2013

diffLL2K = function(image, bandwidth, plot = FALSE){
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
  if (length(bandwidth) != 1)
    stop("bandwidth must be a positive integer")
  if (n1 + 2 * bandwidth + 2 > 600)
    stop("bandwidth is too large or the resolution
of the image is too high.")
  n1 = dim(image)[1]
  z = matrix(as.double(image), ncol = n1)
  k = as.integer(bandwidth)
  out = .Fortran('ll2k_diff', n = as.integer(n1 - 1), obsImg = z,
    bandwidth = as.integer(k), diff = z)
  if (plot == FALSE)
  return(out$diff)
  else x = seq(0, 1, length = n1); y = x
  image(x, y, out$diff, col = gray((0:255)/255))
  return(out$diff)
}
