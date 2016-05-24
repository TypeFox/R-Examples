# This is R source code for function 'modify1' in the R
# package 'image'.
# Creator: Yicheng Kang
# Date: April 30, 2013

modify1 = function(bandwidth, image, edge, plot = FALSE){
  obsImg = image
   if (!is.matrix(obsImg) | dim(obsImg)[1] != dim(obsImg)[2])
    stop('obsImg must be a square matrix') 
  if (!is.matrix(edge) | dim(edge)[1] != dim(edge)[2])
    stop('edge must be a square matrix')
  if (length(edge[(edge != 0) & (edge != 1)]) >= 1)
    stop('edge can only have entry equal to 0 or 1')
  if (dim(obsImg)[1] != dim(edge)[1])
    stop('obsImg and edge must have the same size')
  if (!is.numeric(bandwidth) | length(bandwidth) > 1 |
      as.integer(bandwidth) < 1)
    stop('bandwidth must be a positive integer')
  n1 = dim(edge)[1]
  k1 = as.integer(bandwidth)
  itemp = as.integer(2 * k1 + 1)
  obsImg = matrix(as.double(obsImg), ncol = n1)
  z = array(as.double(0), c(601, 601))
  z[1:n1, 1:n1] = obsImg
  ext = .Fortran('extend', n = as.integer(n1 - 1), k = k1,
    z = z, z1 = array(as.double(0), c(601, 601)))
  edge = matrix(as.integer(edge), ncol = n1)
  edge.ext = array(as.integer(0), c(601, 601))
  edge.ext[k1:(k1 + n1 - 1), k1:(k1 + n1 - 1)] = edge
  out = .Fortran('modify1', n = as.integer(n1 - 1), k =
    itemp, bound = k1, z = ext$z1, edge = edge.ext)  
  edge1 = out$edge[k1:(k1 + n1 - 1), k1 : (k1 + n1 - 1)]
  if (plot == FALSE){
    return(edge1)
  }
  else {
    x = seq(0, 1, length = n1); y = x
    image(x, y, 1 - edge, col = gray(0:1), main = 'Original')
    image(x, y, 1 - edge1, col = gray(0:1), main = 'Modified 1')
    return(edge1)
  }
}
  
