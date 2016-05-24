absstddif <-
function(X_mat, t_ind, std_dif) {
     n_vrs = ncol(X_mat)
     n_obs = nrow(X_mat)
     mom_tols_out = NA
     for (j in 1:n_vrs) {
          yes_before = unlist(X_mat[t_ind == 1, j])
          no_before = unlist(X_mat[t_ind == 0, j])
          pooled_sd = sqrt((var(yes_before, na.rm = TRUE) + var(no_before, na.rm = TRUE))/2)
          mom_tols_out = c(mom_tols_out, pooled_sd*std_dif)
     }
     mom_tols_out = mom_tols_out[-1]
     mom_tols_out
}
