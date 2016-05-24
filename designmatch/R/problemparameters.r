#! Generate the parameters for bpmatch
.problemparameters = function(t_ind, dist_mat, subset_weight, n_controls, total_pairs,
                             mom_covs, mom_tols,
                             ks_covs, ks_n_grid, ks_tols,
                             exact_covs,
                             near_exact_covs, near_exact_devs,
                             fine_covs,
                             near_fine_covs, near_fine_devs,
                             near_covs, near_pairs, near_groups,
                             far_covs, far_pairs, far_groups,
                             use_controls,
                             approximate) {
  
  #! Number of treated units and controls
  n_t = sum(t_ind)
  n_c = length(t_ind)-n_t
  
  #! Number of dec. vars.
  n_dec_vars = n_t*n_c
  
  #! Number of moment covariates
  n_mom_covs = 0
  if(!is.null(mom_covs)) {
    n_mom_covs = ncol(mom_covs)
  }

  #! Number of K-S covariates
  n_ks_covs = 0
  if(!is.null(ks_covs)) {
    n_ks_covs = ncol(ks_covs)
    
    if ((length(ks_n_grid)==1) && (n_ks_covs > 1)) {
      ks_n_grid = rep(ks_n_grid, n_ks_covs)
    }
    
  }
  #! Parameters used to minimize the K-S statistic
  ks_covs_aux = NULL
  if (is.null(ks_covs)) {
    max_ks_n_grid = 0
  }
  if (!is.null(ks_covs)) {
    max_ks_n_grid = max(ks_n_grid)
    #! Grid of values	
    ks_grid = matrix(0, nrow = max_ks_n_grid, ncol = n_ks_covs)
    for (i in 1:n_ks_covs) {
      ks_covs_t_aux = ks_covs[, i][t_ind==1]
      ks_grid_aux = quantile(ks_covs_t_aux, probs = seq(1/ks_n_grid[i], 1, 1/ks_n_grid[i]))
      ks_grid_aux = c(ks_grid_aux, rep(0, max_ks_n_grid-ks_n_grid[i]))
      ks_grid[, i] = ks_grid_aux
    }	
    #! Auxiliary covariates
    ks_covs_aux = matrix(0, nrow = length(t_ind), ncol = max_ks_n_grid*n_ks_covs)
    for (i in 1:n_ks_covs) {
      k = (i-1)*max_ks_n_grid
      for (j in 1:max_ks_n_grid) {
        ks_covs_aux[, j+k][ks_covs[, i]<ks_grid[j, i]] = 1
      }
    }
  }
  
  #! Coeffs. of the obj. fun., cvec
  if (is.null(dist_mat)) {
    if (approximate == 1 | n_controls == 1) {
      cvec = -(1*rep(1, n_t*n_c))
    }
    else {
      cvec = c(-(1*rep(1, n_t*n_c)), rep(0, n_t)) 
    }
  }
  if (!is.null(dist_mat)) {
    if (approximate == 1 | n_controls == 1) {
      cvec = as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c))
    }
    else {
      cvec = c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)), rep(0, n_t))
    }
    
  }
  
  #! Constraint matrix, Amat
  constraintmat_out = .constraintmatrix(t_ind, n_controls, total_pairs,
                                       mom_covs, mom_tols,
                                       ks_covs, ks_covs_aux, ks_n_grid, ks_tols,
                                       exact_covs,
                                       near_exact_covs, near_exact_devs,
                                       fine_covs,
                                       near_fine_covs, near_fine_devs,
                                       near_covs, near_pairs, near_groups,
                                       far_covs, far_pairs, far_groups,
                                       use_controls,
                                       approximate)
  
  cnstrn_mat = constraintmat_out$cnstrn_mat
  bvec_7 = constraintmat_out$bvec_7
  bvec_8 = constraintmat_out$bvec_8
  rows_ind_far_pairs = constraintmat_out$rows_ind_far_pairs
  rows_ind_near_pairs = constraintmat_out$rows_ind_near_pairs
  
  # Constraint vector, bvec
  #! Parts 1 and 2
  if (approximate == 1 | n_controls == 1) {
    bvec = c(rep(n_controls, n_t), rep(1, n_c))
  }
  else {
    bvec = c(rep(0, n_t), rep(1, n_c))
  }
  
  #! Part 3: moments
  if (!is.null(mom_covs)) {
    bvec = c(bvec, rep(0, 2*n_mom_covs))	
  }	
  
  #! Part 4: K-S
  if (!is.null(ks_covs)) {
    bvec = c(bvec, rep(0, 2*n_ks_covs*max_ks_n_grid))
  }
  
  #! Part 5: exact
  if (!is.null(exact_covs)) {
    bvec = c(bvec, rep(0, ncol(exact_covs))) 
  }	
  
  #! Part 6: near-exact
  if (!is.null(near_exact_covs)) {
    bvec = c(bvec, near_exact_devs) 
  }
  
  #! Part 7: fine
  if (!is.null(fine_covs)) {
    bvec = c(bvec, bvec_7) 
  }
  
  #! Part 8: near-fine
  if (!is.null(near_fine_covs)) {
    n_near_fine_covs = ncol(near_fine_covs)
    near_fine_devs_aux = NULL
    for (j in 1:n_near_fine_covs) {
      near_fine_cov = near_fine_covs[, j]
      near_fine_devs_aux = c(near_fine_devs_aux, rep(near_fine_devs[j], length(names(table(near_fine_cov))) ))
    }
    bvec_8_aux = rep(NA, length(bvec_8)*2)
    bvec_8_aux[1:length(bvec_8)] = -near_fine_devs_aux
    bvec_8_aux[(length(bvec_8)+1):(2*length(bvec_8))] = near_fine_devs_aux
    bvec = c(bvec, bvec_8_aux) 
  }
  
  #! Part 9: far
  if (!is.null(far_covs)) {
    n_far_covs = ncol(far_covs)
    for (j in 1:n_far_covs) {
      if (!is.null(far_groups)) {
        bvec = c(bvec, rep(0, 1))
      }
      if (!is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        bvec = c(bvec, rep(0, length(table(rows_ind_far_pairs[[j]]))))
      }
    }
  }
  
  #! Part 10: near
  if (!is.null(near_covs)) {
    n_near_covs = ncol(near_covs)
    for (j in 1:n_near_covs) {
      if (!is.null(near_groups)) {
        bvec = c(bvec, rep(0, 1))
      }
      if (!is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        bvec = c(bvec, rep(0, length(table(rows_ind_near_pairs[[j]]))))
      }
    }
  }
  
  #! Part 11: use controls
  if (!is.null(use_controls)) {
    bvec = c(bvec, sum(use_controls)) 
  }
  
  #! Part 12: total_pairs
  if (!is.null(total_pairs)) {
    bvec = c(bvec, total_pairs)
  }
  
  # Upper bounds, ub
  #! Parts 1 and 2
  #! Part 3: moments
  #! Part 4: K-S
  if (approximate == 1 | n_controls == 1) {
    ub = rep(1, n_t*n_c)
  }	
  else {
    ub = c(rep(1, n_t*n_c), rep(1, n_t))
  }
  
  # Sense, sense
  #! Parts 1 and 2
  #! Part 3: moments
  #! Part 4: K-S
  if (approximate == 1 | n_controls == 1) {
    sense = c(rep("L", n_t), rep("L", n_c), rep("L", 2*n_mom_covs), rep("L", 2*n_ks_covs*max_ks_n_grid))
  }
  else {
    sense = c(rep("E", n_t), rep("L", n_c), rep("L", 2*n_mom_covs), rep("L", 2*n_ks_covs*max_ks_n_grid))
  }
  
  #! Part 5: exact
  if (!is.null(exact_covs)) {
    sense = c(sense, rep("E", ncol(exact_covs))) 
  }
  
  #! Part 6: near-exact
  if (!is.null(near_exact_covs)) {
    sense = c(sense, rep("L", ncol(near_exact_covs))) 
  }
  
  #! Part 7: fine
  if (!is.null(fine_covs)) {
    sense = c(sense, rep("E", length(bvec_7))) 
  }
  
  #! Part 8: near-fine
  if (!is.null(near_fine_covs)) {
    #sense = c(sense, rep(c("G", "L"), length(bvec_8))) 
    sense = c(sense, rep("G", length(bvec_8)), rep("L", length(bvec_8)))
  }
  
  #! Part 9: far
  if (!is.null(far_covs)) {
    n_far_covs = ncol(far_covs)
    for (j in 1:n_far_covs) {
      if (!is.null(far_groups)) { 
        sense = c(sense, rep("G", 1)) 
      }
      if (!is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        sense = c(sense, rep("E", length(table(rows_ind_far_pairs[[j]]))))
      }
    }
  }
  
  #! Part 10: near
  if (!is.null(near_covs)) {
    n_near_covs = ncol(near_covs)
    for (j in 1:n_near_covs) {
      if (!is.null(near_groups)) { 
        sense = c(sense, rep("L", 1)) 
      }
      if (!is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        sense = c(sense, rep("E", length(table(rows_ind_near_pairs[[j]]))))
      }
    }
  }
  
  #! Part 11: use controls
  if (!is.null(use_controls)) {
    sense = c(sense, "E") 
  }
  
  #! Part 12: total_pairs
  if (!is.null(total_pairs)) {
    sense = c(sense, "E")
  }
  
  # Variable types, vtype
  #! Parts 1 and 2
  #! Part 3: moments
  #! Part 4: K-S
  if (approximate == 1) {
    vtype = rep("C", n_t*n_c)
  }
  else if (n_controls == 1) {
    vtype = rep("B", n_t*n_c)
  }
  else {
    vtype = c(rep("B", n_t*n_c), rep("B", n_t))
  }
  
  #! Part 5: exact
  #! Part 6: near-exact
  #! Part 7: fine
  #! Part 8: near-fine
  #! Part 9: far
  #! Part 10: near
  #! Part 11: use controls
  #! Part 12: total_pairs
  
  # c_index
  c_index = rep(1:n_c, n_t)	
  
  # Output
  return(list(n_t = n_t, n_c = n_c,  
              cvec = cvec, 
              Amat = cnstrn_mat, 
              bvec = bvec, 
              ub = ub, 
              sense = sense,
              vtype = vtype,
              c_index = c_index))
}