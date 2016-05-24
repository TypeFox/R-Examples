pstau <- function(matrix, matrix1, num_cat) {
  #Function computes the pseudo tau coefficient (based on the mutability).  
  (nnWCM(matrix1, num_cat) - nnWCM(matrix, num_cat))/nnWCM(matrix1, num_cat)
}