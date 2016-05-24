#'PI_2
#'@keywords internal
PI_1 <- function(t, all_times, gamma_vec, U){

  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time_mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- matrix(rep(1:length(all_times), length(t)), length(all_times), length(t))*(event_time_mat >= time_mat)

  almost_res <- apply(as.matrix(as.vector(exp(as.matrix(U[indexes, ])%*%gamma_vec))*U[indexes, ], ncol = length(gamma_vec)), 2, cumsum)[cumsum(colSums(event_time_mat >= time_mat)), ]
  almost_res <- as.matrix(almost_res, nrow = length(t))
  res <- almost_res - rbind(0, almost_res[-dim(almost_res)[1], , drop = F])
  # discards the result for the -1
  return(matrix(res[-1,], nrow = length(t) - 1))
}