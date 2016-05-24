#'PI_0
#'@keywords internal
PI_0 <- function(t, all_times, gamma_vec, U){

  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time_mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- matrix(rep(1:length(all_times), length(t)), length(all_times), length(t))*(event_time_mat >= time_mat)

  almost_res <- as.vector(cumsum(exp(as.matrix(U[indexes, ])%*%gamma_vec))[cumsum(colSums(event_time_mat >= time_mat))])
  res <- almost_res - c(0, almost_res[-length(almost_res)])
  # discards the result for the -1
  return(res[-1])
}