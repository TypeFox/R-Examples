#'Ahat computation
#'@keywords internal
# The idea behind this calculation is:
# dEN_i(t) is estimated by the mean number of failure times before i.e. multiply by 1/n (typically dN_i(t) = 1 only for 1 out of n times)
# The second derivative of the log partial likelihood
Ahat <- function(all_times, failures, gamma_vec, U){
  # be careful with the definitions here
  n <- length(all_times)

  fail_times <- all_times[failures == 1]

  pi_0 <- PI_0(fail_times, all_times, gamma_vec, U)/n
  pi_1 <- PI_1(fail_times, all_times, gamma_vec, U)/n

  l <- dim(as.matrix(pi_1))[2]
  pi_1_x2 <- as.matrix(pi_1)[,rep(1:l, l)]*as.matrix(pi_1)[,rep(1:l, rep(l,l))]

  pi_2 <- PI_2(fail_times, all_times, gamma_vec, U)/n

  res <- colSums(pi_0^(-2)*(pi_2*pi_0 - pi_1_x2))

  return(matrix(as.vector(t(res)), ncol = ncol(U)))
}

