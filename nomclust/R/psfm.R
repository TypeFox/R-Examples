psfm <- function(matrix, matrix1, num_cat) {
  #Function computes the pseudo F coefficient based on the Mutability
  dim <- dim(matrix)
  #number of variables
  m <- dim[3]
  #number of clusters
  k <- dim[2]
  #number of elements
  n <- sum(matrix[,,1])
  
  nominator <- (n-k)*(nnWCM(matrix1, num_cat)-nnWCM(matrix, num_cat))
  denominator <- (k-1)*nnWCM(matrix, num_cat)
  
  nominator/denominator
}