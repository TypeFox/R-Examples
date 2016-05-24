psu <- function(matrix, matrix1, num_cat) {
  #Function computes the pseudo tau coefficient (based on the entropy).
  (nnWCE(matrix1, num_cat) - nnWCE(matrix, num_cat))/nnWCE(matrix1, num_cat)
}