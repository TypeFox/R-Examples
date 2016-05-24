cpploss <- function(default_distr_a, link_function_a, S_a, Sigma_a, W_a, PD_a, PL_a, calc_rc_a, loss_thr_a, max_entries_a) {
  .Call('GCPM_cpploss', PACKAGE = 'GCPM', default_distr_a, link_function_a, S_a, Sigma_a, W_a, PD_a, PL_a, calc_rc_a, loss_thr_a, max_entries_a)
}
