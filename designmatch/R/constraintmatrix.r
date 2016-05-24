#! Builds the constraint matrix
.constraintmatrix = function(t_ind, n_controls, total_pairs,
                            mom_covs, mom_tols,
                            ks_covs, ks_covs_aux, ks_n_grid, ks_tols,
                            exact_covs,
                            near_exact_covs, near_exact_devs,
                            fine_covs,
                            near_fine_covs, near_fine_devs,
                            near_covs, near_pairs, near_groups,
                            far_covs, far_pairs, far_groups,
                            use_controls,
                            approximate) {
  
  #! Number of treated units, number of controls
  n_t = sum(t_ind)
  n_c = length(t_ind)-n_t	
  
  #! Total number of units
  n_tot = n_t*n_c
  
  #! Build parts of the constraint matrix
  #! Part 1
  if (approximate == 1 | n_controls == 1) {
    row_ind_1 = sort(rep(1:n_t, n_c))
    col_ind_1 = 1:n_tot
    ones_1 = rep(1, n_tot)
  }
  else {
    row_ind_1 = c(sort(rep(1:n_t, n_c)), 1:n_t)
    col_ind_1 = 1:(n_tot+n_t)
    ones_1 = c(rep(1, n_tot), rep(-1*n_controls, n_t))
  }
  
  #! Part 2
  row_ind_2 = sort(rep(1:n_c, n_t))+n_t
  col_ind_2 = rep(seq(1, n_t*n_c, n_c), n_c)+(sort(rep(1:n_c, n_t))-1)
  ones_2 = rep(1, n_tot)
  #! Current max row index
  row_ind_cur	= max(row_ind_2)
  
  #! Parts 3 and 4: moments and K-S
  mom_ks_covs = NULL
  if (!is.null(mom_covs) | !is.null(ks_covs)) {
    row_ind_3.4 = 0
    #! Number of moment covariates
    n_mom_covs = 0
    if(!is.null(mom_covs)) {
      n_mom_covs = ncol(mom_covs)
    }
    #! Number of K-S covariates
    n_ks_covs = 0
    if(!is.null(ks_covs)) {
      n_ks_covs = ncol(ks_covs)
    }
    # Bind moment and K-S covariates		
    if(!is.null(mom_covs) & is.null(ks_covs_aux)) {
      mom_ks_covs = mom_covs
      mom_ks_tols = mom_tols
    }
    if(is.null(mom_covs) & !is.null(ks_covs_aux)) {
      mom_ks_covs = ks_covs_aux
      mom_ks_tols = NA
      for (i in 1:ncol(ks_covs)) {
        mom_ks_tols = c(mom_ks_tols, rep(ks_tols[i], ks_n_grid[i]), rep(0, max(ks_n_grid)-ks_n_grid[i]))
      }
      mom_ks_tols = mom_ks_tols[-1]
    }
    if(!is.null(mom_covs) & !is.null(ks_covs_aux)) {
      mom_ks_covs = cbind(mom_covs, ks_covs_aux)
      mom_ks_tols = mom_tols
      for (i in 1:ncol(ks_covs)) {
        mom_ks_tols = c(mom_ks_tols, rep(ks_tols[i], ks_n_grid[i]), rep(0, max(ks_n_grid)-ks_n_grid[i]))
      }
    }	
  }					 
  if (!is.null(mom_ks_covs)) {
    n_mom_ks_covs = ncol(mom_ks_covs)
    if (!is.null(mom_tols) | !is.null(ks_tols)) {
      row_ind_3.4 = sort(rep(1:(2*n_mom_ks_covs)+n_t+n_c, n_tot))
    }	
    col_ind_3.4 = NA
    mom_ks_vals_3.4 = NA
    j = 1
    k = 0
    for (i in 1:n_mom_ks_covs) {
      if (n_mom_covs != 0 & i <= n_mom_covs) {
        if (!is.null(mom_tols) | !is.null(ks_tols)) {
          col_ind_3.4 = c(col_ind_3.4, rep(1:n_tot, 2))
        }
      }
      if (n_ks_covs != 0 & i > n_mom_covs) {
        if (!is.null(mom_tols) | !is.null(ks_tols)) {
          col_ind_3.4 = c(col_ind_3.4, rep(1:n_tot, 2))
          k = k+1
          if (k >= max(ks_n_grid)) {
            j = j+1
            k = 0	
          }
        }
      }
      temp_mean_1 = rep(mom_ks_covs[t_ind==0, i], n_t)-(mom_ks_covs[t_ind==1, i])[sort(rep(1:n_t, n_c))]
      if (!is.null(mom_tols) | !is.null(ks_tols)) {
        temp_mean_2 = temp_mean_1-(mom_ks_tols[i]*rep(1, n_t*n_c))
        temp_mean_3 = -temp_mean_1-(mom_ks_tols[i]*rep(1, n_t*n_c))
      }	
      mom_ks_vals_3.4 = c(mom_ks_vals_3.4, temp_mean_2, temp_mean_3)
      if (i == 1) {
        col_ind_3.4 = col_ind_3.4[-1]
        mom_ks_vals_3.4 = mom_ks_vals_3.4[-1]
      }
    }
    #! Current max row index
    row_ind_cur	= max(row_ind_3.4)
  }	
  
  #! Part 5: exact
  rows_exact = NULL
  cols_exact = NULL
  vals_exact = NULL
  if (!is.null(exact_covs)) {
    n_exact_cats = ncol(exact_covs)
    for (i in 1:n_exact_cats) {
      rows_exact = c(rows_exact, rep(row_ind_cur+i, n_t*n_c))
      cols_exact = c(cols_exact, 1:(n_t*n_c))
      dist_exact_cov = abs(outer(exact_covs[t_ind==1, i], exact_covs[t_ind==0, i], "-"))
      dist_exact_cov = t(dist_exact_cov)
      vals_exact = c(vals_exact, as.vector(dist_exact_cov))
    }	
    row_ind_5 = rows_exact
    col_ind_5 = cols_exact
    exact_vals_5 = vals_exact
    row_ind_cur	= max(row_ind_5)
  }
  
  #! Part 6: near-exact
  rows_near_exact = NULL
  cols_near_exact = NULL
  vals_near_exact = NULL
  if (!is.null(near_exact_covs)) {
    n_near_exact_cats = ncol(near_exact_covs)
    for (i in 1:n_near_exact_cats) {
      rows_near_exact = c(rows_near_exact, rep(row_ind_cur+i, n_t*n_c))
      cols_near_exact = c(cols_near_exact, 1:(n_t*n_c))
      dist_near_exact_cov = abs(outer(near_exact_covs[t_ind==1, i], near_exact_covs[t_ind==0, i], "-"))
      dist_near_exact_cov = t(dist_near_exact_cov)
      vals_near_exact = c(vals_near_exact, as.vector(dist_near_exact_cov))
    }	
    row_ind_6 = rows_near_exact
    col_ind_6 = cols_near_exact
    near_exact_vals_6 = vals_near_exact
    row_ind_cur	= max(row_ind_6)
  }
  
  #! Part 7: fine
  bvec_7 = NULL
  rows_fine = NULL
  cols_fine = NULL
  vals_fine = NULL
  if (!is.null(fine_covs)) {
    #! Transform fine_covs to a matrix of binary inds. for each cat. of each fine balancing covariate
    fine_covs_2 = rep(NA, nrow(fine_covs))
    n_fine_covs = ncol(fine_covs)
    j = 1
    for (i in 1:n_fine_covs) {	
      aux = factor(fine_covs[, i])
      fine_covs_2 = cbind(fine_covs_2, diag(nlevels(aux))[aux,])
      if (j == 1) {
        fine_covs_2 = fine_covs_2[, -1]
      }
      j = j+1
    }
    n_fine_cats = ncol(fine_covs_2)
    for (i in 1:n_fine_cats) {
      rows_fine = c(rows_fine, rep(row_ind_cur+i, n_t*n_c))
      cols_fine = c(cols_fine, 1:(n_t*n_c))
      dist_fine_cov = outer(fine_covs_2[t_ind==1, i], fine_covs_2[t_ind==0, i], "-")
      dist_fine_cov = t(dist_fine_cov)
      vals_fine = c(vals_fine, as.vector(dist_fine_cov))
    }	
    row_ind_7 = rows_fine
    col_ind_7 = cols_fine
    fine_vals_7 = vals_fine
    bvec_7 = rep(0, n_fine_cats)
    row_ind_cur	= max(row_ind_7)
  }
  
  #! Part 8: near-fine
  bvec_8 = NULL
  rows_near_fine = NULL
  cols_near_fine = NULL
  vals_near_fine = NULL
  if (!is.null(near_fine_covs)) {
    #! Transform fine_covs to a matrix of binary inds. for each cat. of each fine balancing covariate
    near_fine_covs_2 = rep(NA, nrow(near_fine_covs))
    n_near_fine_covs = ncol(near_fine_covs)
    j = 1
    for (i in 1:n_near_fine_covs) {  
      near_aux = factor(near_fine_covs[, i])
      near_fine_covs_2 = cbind(near_fine_covs_2, diag(nlevels(near_aux))[near_aux,])
      if (j == 1) {
        near_fine_covs_2 = near_fine_covs_2[, -1]
      }
      j = j+1
    }
    n_near_fine_cats = ncol(near_fine_covs_2)
    for (i in 1:n_near_fine_cats) {
      rows_near_fine = c(rows_near_fine, rep(row_ind_cur+i, n_t*n_c))
      cols_near_fine = c(cols_near_fine, 1:(n_t*n_c))
      dist_near_fine_cov = outer(near_fine_covs_2[t_ind==1, i], near_fine_covs_2[t_ind==0, i], "-")
      dist_near_fine_cov = t(dist_near_fine_cov)
      vals_near_fine = c(vals_near_fine, as.vector(dist_near_fine_cov))
    }
    row_ind_cur = max(rows_near_fine)
    for (i in 1:n_near_fine_cats) {
      rows_near_fine = c(rows_near_fine, rep(row_ind_cur+i, n_t*n_c))
    }
    row_ind_8 = rows_near_fine
    col_ind_8 = c(cols_near_fine, cols_near_fine)
    near_fine_vals_8 = c(vals_near_fine, vals_near_fine)
    bvec_8 = rep(0, n_near_fine_cats)
    row_ind_cur	= max(row_ind_8)
  }
  
  #! Part 9: Far
  rows_ind_far_pairs = list()
  if (!is.null(far_covs)) {
    row_ind_9 = NULL
    col_ind_9 = NULL
    far_cov_vals_9 = NULL
    n_far_covs = ncol(far_covs)
    for (j in 1:n_far_covs) {
      far_cov = far_covs[,j]
      #! Far on average constraints
      if (!is.null(far_groups)) {
        far_group = far_groups[j]
        row_ind_far_all = sort(c(rep(row_ind_cur+1, n_tot)))		
        col_ind_far_all = rep(1:n_tot, 1)			
        temp_mean_3 = (-rep(far_cov[t_ind==0], n_t)+((far_cov[t_ind==1])[sort(rep(1:n_t, n_c))]))-(far_group*rep(1, n_t*n_c))
        vals_far_all = c(temp_mean_3)		
        row_ind_cur	= max(row_ind_far_all)
      }
      #! Far on all pairs constraints
      if (!is.null(far_pairs)) {
        far_pair = far_pairs[j]
        aux = abs(outer(far_cov[t_ind==1], far_cov[t_ind==0], FUN = "-"))
        temp = as.vector(matrix(t(aux), nrow = 1, byrow = TRUE))
        cols_ind_far_pairs = which(temp<far_pair)
        if (length(cols_ind_far_pairs)>0) {
          rows_ind_far_pairs[[j]] = row_ind_cur+(1:length(cols_ind_far_pairs))
          vals_far_pairs = rep(1, length(cols_ind_far_pairs))
          row_ind_cur	= max(rows_ind_far_pairs[[j]])
        }
        if (length(cols_ind_far_pairs)==0) {
          cols_ind_far_pairs = NULL
          rows_ind_far_pairs[[j]] = -1
          vals_far_pairs = NULL
        }
      }		
      #! Put together
      if (!is.null(far_groups) && is.null(far_pairs)) {
        row_ind_9 = c(row_ind_9, row_ind_far_all)
        col_ind_9 = c(col_ind_9, col_ind_far_all)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_all)
      }
      if (is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        row_ind_9 = c(row_ind_9, rows_ind_far_pairs[[j]])
        col_ind_9 = c(col_ind_9, cols_ind_far_pairs)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_pairs)
      }
      if (!is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        row_ind_9 = c(row_ind_9, row_ind_far_all, rows_ind_far_pairs[[j]])
        col_ind_9 = c(col_ind_9, col_ind_far_all, cols_ind_far_pairs)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_all, vals_far_pairs)
      }
      if (!is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] == -1) {
        row_ind_9 = c(row_ind_9, row_ind_far_all)
        col_ind_9 = c(col_ind_9, col_ind_far_all)
        far_cov_vals_9 = c(far_cov_vals_9, vals_far_all)
      }
    }
  }
  
  #! Part 10: Near
  rows_ind_near_pairs = list()
  if (!is.null(near_covs)) {
    row_ind_10 = NULL
    col_ind_10 = NULL
    near_cov_vals_10 = NULL
    n_near_covs = ncol(near_covs)
    for (j in 1:n_near_covs) {
      near_cov = near_covs[,j]
      #! Near on average constraints
      if (!is.null(near_groups)) {
        near_group = near_groups[j]
        row_ind_near_all = sort(c(rep(row_ind_cur+1, n_tot)))  	
        col_ind_near_all = rep(1:n_tot, 1)			
        temp_mean_4 = (-rep(near_cov[t_ind==0], n_t)+((near_cov[t_ind==1])[sort(rep(1:n_t, n_c))]))-(near_group*rep(1, n_t*n_c))
        vals_near_all = c(temp_mean_4)		
        row_ind_cur	= max(row_ind_near_all)
      }
      #! Near on all pairs constraints
      if (!is.null(near_pairs)) {
        near_pair = near_pairs[j]
        aux = abs(outer(near_cov[t_ind==1], near_cov[t_ind==0], FUN = "-"))
        temp = as.vector(matrix(t(aux), nrow = 1, byrow = TRUE))
        cols_ind_near_pairs = which(temp>near_pair)
        if (length(cols_ind_near_pairs)>0) {
          rows_ind_near_pairs[[j]] = row_ind_cur+(1:length(cols_ind_near_pairs))
          vals_near_pairs = rep(1, length(cols_ind_near_pairs))
          row_ind_cur	= max(rows_ind_near_pairs[[j]])
        }
        if (length(cols_ind_near_pairs)==0) {
          cols_ind_near_pairs = NULL
          rows_ind_near_pairs[[j]] = -1
          vals_near_pairs = NULL
        }
      }		
      #! Put together
      if (!is.null(near_groups) && is.null(near_pairs)) {
        row_ind_10 = c(row_ind_10, row_ind_near_all)
        col_ind_10 = c(col_ind_10, col_ind_near_all)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_all)
      }
      if (is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        row_ind_10 = c(row_ind_10, rows_ind_near_pairs[[j]])
        col_ind_10 = c(col_ind_10, cols_ind_near_pairs)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_pairs)
      }
      if (!is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        row_ind_10 = c(row_ind_10, row_ind_near_all, rows_ind_near_pairs[[j]])
        col_ind_10 = c(col_ind_10, col_ind_near_all, cols_ind_near_pairs)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_all, vals_near_pairs)
      }
      if (!is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] == -1) {
        row_ind_10 = c(row_ind_10, row_ind_near_all)
        col_ind_10 = c(col_ind_10, col_ind_near_all)
        near_cov_vals_10 = c(near_cov_vals_10, vals_near_all)
      }
    }
  }
  
  # Part 11: use controls
  if (!is.null(use_controls)) {
    use_controls = use_controls[(n_t+1):(n_t+n_c)]
    use_controls_aux = rep(use_controls, n_t)
    col_ind_11 = (1:n_tot)[use_controls_aux==1]
    row_ind_11 = rep(row_ind_cur+1, length(col_ind_11))
    use_controls_vals_11 = rep(1, length(col_ind_11))
    
    row_ind_cur = max(row_ind_11)
  }
  
  # Part 12: total_pairs
  if (!is.null(total_pairs)) {
    row_ind_12 = rep(row_ind_cur+1, n_t*n_c)
    col_ind_12 = 1:(n_t*n_c)
    ones_12 = rep(1, n_t*n_c)
    
    row_ind_cur = max(row_ind_12)
  }
  
  #! Put all the parts of the constraint matrix together
  #! Parts 1 and 2
  row_ind = c(row_ind_1, row_ind_2)
  col_ind = c(col_ind_1, col_ind_2)
  vals = c(ones_1, ones_2)
  #! Parts 3 and 4
  if (!is.null(mom_ks_covs)) {
    row_ind = c(row_ind, row_ind_3.4)
    col_ind = c(col_ind, col_ind_3.4)
    vals = c(vals, mom_ks_vals_3.4)
  }
  #! Part 5
  if (!is.null(exact_covs)) {
    row_ind = c(row_ind, row_ind_5)
    col_ind = c(col_ind, col_ind_5)
    vals = c(vals, exact_vals_5)
  }
  #! Part 6
  if (!is.null(near_exact_covs)) {
    row_ind = c(row_ind, row_ind_6)
    col_ind = c(col_ind, col_ind_6)
    vals = c(vals, near_exact_vals_6)
  }
  #! Part 7
  if (!is.null(fine_covs)) {
    row_ind = c(row_ind, row_ind_7)
    col_ind = c(col_ind, col_ind_7)
    vals = c(vals, fine_vals_7)
  }
  #! Part 8
  if (!is.null(near_fine_covs)) {
    row_ind = c(row_ind, row_ind_8)
    col_ind = c(col_ind, col_ind_8)
    vals = c(vals, near_fine_vals_8)
  }
  #! Part 9
  if (!is.null(far_covs)) {
    row_ind = c(row_ind, row_ind_9)
    col_ind = c(col_ind, col_ind_9)
    vals = c(vals, far_cov_vals_9)	
  }	
  #! Part 10
  if (!is.null(near_covs)) {
    row_ind = c(row_ind, row_ind_10)
    col_ind = c(col_ind, col_ind_10)
    vals = c(vals, near_cov_vals_10)  
  }	
  #! Part 11				 
  if (!is.null(use_controls)) {
    row_ind = c(row_ind, row_ind_11)
    col_ind = c(col_ind, col_ind_11)
    vals = c(vals, use_controls_vals_11)	
  }
  #! Part 12
  if (!is.null(total_pairs)) {
    row_ind = c(row_ind, row_ind_12)
    col_ind = c(col_ind, col_ind_12)
    vals = c(vals, ones_12)
  }
  
  aux = cbind(row_ind, col_ind, vals)[order(col_ind), ]
  cnstrn_mat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  #! Output  
  return(list(cnstrn_mat = cnstrn_mat, bvec_7 = bvec_7, bvec_8 = bvec_8, rows_ind_far_pairs = rows_ind_far_pairs, rows_ind_near_pairs = rows_ind_near_pairs))
    
}