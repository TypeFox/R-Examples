#'Lambda and its perturbation
#'
#'@keywords internal
#'@rdname lambda
lambda <- function(t, all_times, failures, gamma_vec, U){

  fail_times <- all_times[failures == 1]

  sums_fail_times <- PI_0(fail_times, all_times, gamma_vec, U)
  denoms <- 1/sums_fail_times

  time_mat <- matrix(rep(t, length(fail_times)), length(fail_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(fail_times, length(t)), length(fail_times), length(t))

  return(t(event_time_mat <= time_mat)%*%denoms)
}


#'@keywords internal
#'@rdname lambda
lambda_pert <- function(t, perturb_mat, all_times, failures, gamma_vec, U){
  fail <- which(failures == 1)
  fail_times <- all_times[fail]

  time_mat <- matrix(rep(t, length(fail_times)), length(fail_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(fail_times, length(t)), length(fail_times), length(t))

  fail_times_less_true_false <- (event_time_mat <= time_mat)

  n <- length(all_times)

  dM_mat <- dM(all_times, failures, gamma_vec, U)
  precomp_exp_sum <- PI_0(fail_times, all_times, gamma_vec, U)
  dM_over_precomp_exp_sum <- dM_mat/precomp_exp_sum

  perturbated_mat_across_fail_times <- t(perturb_mat)%*%t(dM_over_precomp_exp_sum)

  fail_times_less_true_false <- (event_time_mat <= time_mat)

  # this is exactly lambda
  denom_sum_vec <- t(fail_times_less_true_false)%*%(1/precomp_exp_sum)

  perturbated_with_only_needed_obs <- t(fail_times_less_true_false)%*%t(perturbated_mat_across_fail_times)

  # just repeat the lambda vector
  first_part <- matrix(rep(denom_sum_vec, ncol(perturb_mat)), length(t))

  # the sum of SUM(xi_i dM_i(t)/PI_0(t))
  second_part <- perturbated_with_only_needed_obs

  # COMPUTATION OF THE THIRD PART
  res <- exp(log(first_part) + (second_part)/(as.vector(denom_sum_vec)))

  return(res)
}


