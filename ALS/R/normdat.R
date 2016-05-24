"normdat" <-
function (mat) 
{
  if(!is.matrix(mat))
    mat <- as.matrix(mat)
  for (i in 1:ncol(mat))
    mat[, i] <- mat[, i]/max(abs(mat[, i]))
  mat

}

