#'dM
#'@keywords internal
dM <- function(all_times, failures, gamma_vec, U){

  fail_times <- all_times[which(failures == 1)]

  almost_ans <- M_vec(fail_times, all_times, failures, gamma_vec, U)

  ans <- almost_ans - rbind(rep(0, length(all_times)), almost_ans[-nrow(almost_ans), ])

  return(ans)
}