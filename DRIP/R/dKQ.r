# This is R source code for function 'dKQ' in the R package
# 'image'.

dKQ = function(edge1, edge2){
  if (!is.matrix(edge1) | !is.matrix(edge2))
    stop('Both arguments must be matrices.')
  if ((dim(edge1)[1] != dim(edge2)[1]) |
      (dim(edge1)[2] != dim(edge2)[2]))
    stop('Two matrices must be of the same size.')
  if (dim(edge1)[1] != dim(edge1)[2])
    stop('Argument must be a square matrix.')
  if ((length(edge1[(edge1 != 0) & (edge1 != 1)]) >= 1) |
      (length(edge2[(edge2 != 0) & (edge2 != 1)]) >= 1))
    stop('Both matrices can only have entries of 0 or 1.')
  n = dim(edge1)[1]
  out = .Fortran('d_KQ', n = as.integer(n - 1), edge1 =
    matrix(as.integer(edge1), ncol = n), edge2 =
    matrix(as.integer(edge2), ncol = n), dKQ = as.double(100))
  return(out$dKQ)
}

