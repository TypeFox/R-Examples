# This is R source code for function 'modify2' in the R
# package 'image'.
# Creator: Yicheng Kang
# Date: April 29, 2013

modify2 = function(bandwidth, edge, plot = FALSE){
  if (!is.matrix(edge) | dim(edge)[1] != dim(edge)[2])
    stop('edge must be a square matrix')
  if (length(edge[(edge != 0) & (edge != 1)]) >= 1)
    stop('edge can only have entry equal to 0 or 1')
  if (!is.numeric(bandwidth) | length(bandwidth) > 1 |
      as.integer(bandwidth) < 1)
    stop('bandwidth must be a positive integer')
  n1 = dim(edge)[1]
  k1 = as.integer(bandwidth)
  itemp = as.integer(2 * k1 + 1)
  edge = matrix(as.integer(edge), ncol = n1)
  edge.ext = array(as.integer(0), c(601, 601))
  edge.ext[k1:(k1 + n1 - 1), k1:(k1 + n1 - 1)] = edge
  out = .Fortran('modify2', n = as.integer(n1 - 1), k = itemp,
    bound = k1, edge = edge.ext)
  edge2 = out$edge[k1:(k1 + n1 - 1), k1 : (k1 + n1 - 1)]
  if (plot == FALSE){
    return(edge2)
  }
  else {
    x = seq(0, 1, length = n1); y = x
    image(x, y, 1 - edge, col = gray(0:1), main = 'Original')
    image(x, y, 1 - edge2, col = gray(0:1), main = 'Modified 2')
    return(edge2)
  }
}
  
