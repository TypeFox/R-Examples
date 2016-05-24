#'M vector and its perturbation
#'@rdname M_vec
#'@keywords internal
M_vec <- function(t, all_times, failures, gamma_vec, U){

  # the failure times
  fail_times <- all_times[failures == 1]

  # the trick with the 0 infront is done if we have a failure before events
  lambdas_of_fail_times <- c(0, lambda(fail_times, all_times, failures, gamma_vec, U))

  time_mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- (event_time_mat <= time_mat) & (failures == 1)

  # this returns the max index of the max failure time prior to min(t, all_times[i])
  fail_indexes <- try(apply(indexes, 2, cumsum))
  if(inherits(fail_indexes, "try-error")){
    stop("no failures detected: estimation in internal function 'M_vec' impossible")
  }

  # the + 1 here of the indexes is added if we ave failures before events refer to the logic above
  s <- as.vector(exp(U%*%gamma_vec))*(lambdas_of_fail_times[(fail_indexes + 1)])

  ans <- matrix(as.vector(indexes) - s, length(t), length(all_times), byrow = T)
  return(ans)
}


#'M vector perturbation
#'@rdname M_vec
#'@keywords internal
M_vec_pert <- function(perturb_mat, t, all_times, failures, gamma_vec, U){

  # the failure times
  fail_times <- all_times[failures == 1]

  ## comment the last two arguments for log lambda pert
  lambdas_of_fail_times <- rbind(rep(0, ncol(perturb_mat)),
                                 lambda_pert(t=fail_times, perturb_mat, all_times, failures, gamma_vec, U)
                                 )

  time_mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- (event_time_mat <= time_mat) & (failures == 1)

  fail_indexes <- apply(indexes, 2, cumsum)
  if(inherits(fail_indexes, "try-error")){
    stop("no failures detected: estimation in internal function 'M_vec_pert' impossible")
  }

  # the + 1 here of the indexes is added if we have failures before events refer to the logic above
  lambdas_of_fail_times_prev <- c(0, lambda(fail_times, all_times, failures, gamma_vec, U))
  prev_s <- as.vector(exp(U%*%gamma_vec))*(lambdas_of_fail_times_prev[(fail_indexes + 1)])


  gamma_star <- matrix(rep(gamma_vec, ncol(perturb_mat)), length(gamma_vec)) #Compute the original gammas

  exp_matr <- as.matrix(exp(U%*%gamma_star))
  s <- exp_matr[rep(1:nrow(exp_matr), length(t)),]*(lambdas_of_fail_times[(fail_indexes + 1), ])

  ans <- (as.vector(indexes) - s)*perturb_mat[rep(1:nrow(perturb_mat), length(t)),]
  perturbed_part_of_M <- s

  M_v <- matrix(as.vector(indexes) - prev_s, length(t), length(all_times), byrow = T)

  M_v_pert <- as.vector(t(M_v))*(perturb_mat)[rep(1:nrow(perturb_mat), dim(M_v)[1]), ]

  ans <- M_v_pert - (perturbed_part_of_M - prev_s)
  return(ans)
}

