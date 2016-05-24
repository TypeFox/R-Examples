# This is R source code for function 'roofEdge', in the R
# package "image".
# Date: May 6, 2013
# Creator: Yicheng Kang

roofEdge = function(image, bandwidth, thresh, edge1, blur = FALSE,
  plot = FALSE){
  if (!is.matrix(image))
    stop("image data must be a matrix.")
  n1 = as.integer(dim(image)[1])
  n2 = as.integer(dim(image)[2])
  if (n1 != n2)
    stop("image data must be a square matrix.")
  if (!is.numeric(bandwidth) | length(bandwidth) > 1)
    stop("bandwidth must be a positive integer.")
  if (as.integer(bandwidth) < 1)
    stop("bandwidth must be a positive integer.")
  if (n1 + 2 * bandwidth + 2 > 600)
    stop("bandwidth is too large or the resolution
of the image is too high.")
  if (!is.numeric(thresh) | length(thresh) > 1 | thresh < 0)
      stop("threshold  must be a positive number.")
  if(!is.matrix(edge1) | ncol(edge1) != nrow(edge1))
        stop("edge1 must be a square matrix")
  if(!all(edge1==0 | edge1==1))
      stop("edge1's must be either 0 or 1.")
  if(ncol(edge1) != n1)
      stop("edge1 and image are not of the same size.")
  n1 = dim(image)[1]
  z = matrix(as.double(image), ncol = n1)
  edge1 = matrix(as.integer(edge1), ncol=n1)
  edge2 = array(as.integer(0), c(n1, n1))
  k = as.integer(bandwidth)
  u = as.double(thresh)
  if (blur == FALSE) {
    out = .Fortran('roofDetect_denoise', n = as.integer(n1 - 1),
      obsImg = z, bandwidth = as.integer(k), thresh=u, edge1=edge1, edge2=edge2)
  }
  else {
    out = .Fortran('roofDetect_deblur', n = as.integer(n1 - 1),
      obsImg = z, bandwidth = as.integer(k), thresh=u, edge1=edge1, edge2=edge2)
  }
  edge = out$edge2
  if (plot == FALSE) { return(edge) }
  if (plot == TRUE) {
    x = seq(0, 1, length = n1); y = x
    image(x, y, 1 - edge, col = gray(0:1))
    return(edge)
  }
}
