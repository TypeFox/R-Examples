psfe <- function(matrix, matrix1, num_cat) {
  #Function computes the pseudo F coefficient based on the entropy
  dim <- dim(matrix)
  #number of variables
  m <- dim[3]
  #number of clusters
  k <- dim[2]
  #number of elements
  n <- sum(matrix[,,1])

  nominator <- (n-k)*(nnWCE(matrix1, num_cat)-nnWCE(matrix, num_cat))
  denominator <- (k-1)*nnWCE(matrix, num_cat)
  
  nominator/denominator
}
