#'PI_2
#'@keywords internal
PI_2 <- function(t, all_times, gamma_vec, U){
  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time_mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- matrix(rep(1:length(all_times), length(t)), length(all_times), length(t))*(event_time_mat >= time_mat)

  # this is a matrix dim(U)[2]^2 x sum(event_time_mat >= time_mat)
  # it contains the matrices of multiplying the transposed rows of U by themselves (cols of t(U))
  # each column of this matrix is the combined rows of a matrices described above
  l <- ncol(as.matrix(U))
  simplified_mat <- matrix(U[indexes, rep(1:l, l)]*U[indexes, rep(1:l, rep(l,l))], ncol = l^2)

  almost_res <- apply(as.vector(exp(as.matrix(U[indexes, ])%*%gamma_vec))*(simplified_mat), 2, cumsum)[cumsum(colSums(event_time_mat >= time_mat)), ]
  almost_res <- as.matrix(almost_res, ncol = l^2)
  res <- almost_res - rbind(0, as.matrix(almost_res[-dim(almost_res)[1], , drop = F]))
  # discards the result for the -1
  return(matrix(res[-1,], ncol = l^2))
}
