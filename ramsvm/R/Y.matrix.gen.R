Y.matrix.gen <- function(k, kd, nobs, y.train) {

  XI <- XI.gen(k = k, kd = kd)

  Y.matrix <- matrix(data = 0.0, nrow = nobs, ncol = k-1L)

  for( ii in 1L:nobs ) Y.matrix[ii,] <- XI[,y.train[ii]]

  return( Y.matrix )

}
