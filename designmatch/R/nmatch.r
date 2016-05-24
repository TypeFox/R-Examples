nmatch = function(dist_mat, subset_weight = NULL, total_pairs = NULL,
                    mom = NULL,
                    exact = NULL,
                    near_exact = NULL,
                    fine = NULL,
                    near_fine = NULL,
                    near = NULL,
                    far = NULL,
                    solver = NULL) {
  
  if (is.null(mom)) {
    mom_covs = NULL
    mom_tols = NULL
  } else {
    mom_covs = mom$covs
    mom_tols = mom$tols
  }
  
  if (is.null(exact)) {
    exact_covs = NULL
  } else {
    exact_covs = exact$covs
  }
  
  if (is.null(near_exact)) {
    near_exact_covs = NULL
    near_exact_devs = NULL
  } else {
    near_exact_covs = near_exact$covs
    near_exact_devs = near_exact$devs
  }
  
  if (is.null(fine)) {
    fine_covs = NULL
  } else {
    fine_covs = fine$covs
  }
  
  if (is.null(near_fine)) {
    near_fine_covs = NULL
    near_fine_devs = NULL
  } else {
    near_fine_covs = near_fine$covs
    near_fine_devs = near_fine$devs
  }
  
  if (is.null(near)) {
    near_covs = NULL
    near_pairs = NULL
    near_groups = NULL
  } else {
    near_covs = near$covs
    near_pairs = near$pairs
    near_groups = near$groups
  }
  
  if (is.null(far)) {
    far_covs = NULL
    far_pairs = NULL
    far_groups = NULL
  } else {
    far_covs = far$covs
    far_pairs = far$pairs
    far_groups = far$groups
  }
  
  if (is.null(solver)) {
    solver = 'glpk'
    t_max = 60 * 15
    approximate = 1
  } else {
    t_max = solver$t_max
    approximate = solver$approximate
    trace = solver$trace
    round_cplex = solver$round_cplex
    solver = solver$name
  }
  
  cat(format("  Building the matching problem..."), "\n")
  
  #! Total number of units
  n_tot = nrow(dist_mat)
  
  #! Total number of decision variables
  n_dec = (n_tot*(n_tot-1))-sum(1:(n_tot-1))
  
  if (is.null(subset_weight)) {
    subset_weight = 0
  }
  
  #! cvec  
  cvec = t(dist_mat)[lower.tri(dist_mat)]-(subset_weight*rep(1, n_dec))
  
  #! Amat		
  #rows_far_all = NULL
  #cols_far_all = NULL
  #vals_far_all = NULL
  #rows_far_pairs = NULL
  #cols_far_pairs = NULL
  #vals_far_pairs = NULL
  #rows_near_all = NULL
  #cols_near_all = NULL
  #vals_near_all = NULL
  #rows_near_pairs = NULL
  #cols_near_pairs = NULL
  #vals_near_pairs = NULL
  rows_far=  NULL
  cols_far = NULL
  vals_far = NULL
  rows_near = NULL
  cols_near = NULL
  vals_near = NULL
  rows_mom = NULL
  cols_mom = NULL
  vals_mom = NULL
  rows_exact = NULL
  cols_exact = NULL
  vals_exact = NULL
  rows_near_exact = NULL
  cols_near_exact = NULL
  vals_near_exact = NULL
  rows_fine = NULL
  cols_fine = NULL
  vals_fine = NULL
  rows_near_fine = NULL
  cols_near_fine = NULL
  vals_near_fine = NULL
  rows_n = NULL
  cols_n = NULL
  vals_n = NULL		
  
  #! Nonbipartite matching constraints
  rows_nbm = sort(rep(1:n_tot, n_tot-1))
  temp = matrix(0, nrow = n_tot, ncol = n_tot)
  temp[lower.tri(temp)] = 1:n_dec
  temp = temp+t(temp)
  diag(temp) = NA
  cols_nbm = as.vector(t(temp))
  cols_nbm = cols_nbm[!is.na(cols_nbm)]
  vals_nbm = rep(1, (n_tot-1)*n_tot)
  row_count = max(rows_nbm)
  
  #! Far constraints
  rows_ind_far_pairs = list()
  if (!is.null(far_covs)) {
    rows_far = NULL
    cols_far = NULL
    vals_far = NULL
    n_far_covs = ncol(far_covs)
    for (j in 1:n_far_covs) {
      far_cov = far_covs[,j]
      #! Far on average constraints
      if (!is.null(far_groups)) {
        far_group = far_groups[j]
        row_ind_far_all = rep(row_count+1, n_dec)  	
        col_ind_far_all = rep(1:n_dec, 1)			
        i_ind = rep(1:(n_tot-1), (n_tot-1):1)
        aux = matrix(rep(1:n_tot, n_tot), nrow = n_tot, byrow = F)
        j_ind = aux[lower.tri(aux)]
        vals_far_all = 	far_cov[i_ind]-far_cov[j_ind]-(far_group*rep(1, n_dec))	
        row_count	= max(row_ind_far_all)
      }
      #! Far on all pairs constraints
      if (!is.null(far_pairs)) {
        far_pair = far_pairs[j]
        aux = abs(outer(far_cov, far_cov, FUN = "-"))
        temp = as.vector(matrix(t(aux)[lower.tri(aux)], nrow = 1, byrow = TRUE))
        cols_ind_far_pairs = which(temp<far_pair)
        if (length(cols_ind_far_pairs)>0) {
          rows_ind_far_pairs[[j]] = row_count+(1:length(cols_ind_far_pairs))
          vals_far_pairs = rep(1, length(cols_ind_far_pairs))
          row_count	= max(rows_ind_far_pairs[[j]])
        }
        if (length(cols_ind_far_pairs)==0) {
          cols_ind_far_pairs = NULL
          rows_ind_far_pairs[[j]] = -1
          vals_far_pairs = NULL
        }
      }		
      #! Put together
      if (!is.null(far_groups) && is.null(far_pairs)) {
        rows_far = c(rows_far, row_ind_far_all)
        cols_far = c(cols_far, col_ind_far_all)
        vals_far = c(vals_far, vals_far_all)
      }
      if (is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        rows_far = c(rows_far, rows_ind_far_pairs[[j]])
        cols_far = c(cols_far, cols_ind_far_pairs)
        vals_far = c(vals_far, vals_far_pairs)
      }
      if (!is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] != -1) {
        rows_far = c(rows_far, row_ind_far_all, rows_ind_far_pairs[[j]])
        cols_far = c(cols_far, col_ind_far_all, cols_ind_far_pairs)
        vals_far = c(vals_far, vals_far_all, vals_far_pairs)
      }
      if (!is.null(far_groups) && !is.null(far_pairs) && rows_ind_far_pairs[[j]] == -1) {
        rows_far = c(rows_far, row_ind_far_all)
        cols_far = c(cols_far, col_ind_far_all)
        vals_far = c(vals_far, vals_far_all)
      }
    }
  }
  
  #! Far on average constraints
  #if (!is.null(far_covs)) {
  #  rows_far_all = rep(row_count+1, n_dec)
  #  cols_far_all = 1:n_dec
  #  i_ind = rep(1:(n_tot-1), (n_tot-1):1)
  #  aux = matrix(rep(1:n_tot, n_tot), nrow = n_tot, byrow = F)
  #  j_ind = aux[lower.tri(aux)]
  #  vals_far_all = far_covs[i_ind]-far_covs[j_ind]-(far_groups*rep(1, n_dec))
  #  row_count = max(rows_far_all)
  #}
  
  #! Far on all pairs constraints
  #if (!is.null(far_covs)) {
  #  aux = abs(outer(far_covs, far_covs, FUN = "-"))
  #  temp = aux[lower.tri(aux)]
  #  cols_far_pairs = which(temp<far_pairs)
  #  rows_far_pairs = row_count+(1:length(cols_far_pairs))
  #  vals_far_pairs = rep(1, length(cols_far_pairs))
  #  row_count = max(rows_far_pairs)
  #}
  #! Far on all pairs constraints
  #if (!is.null(far_covs)) {
  #  for (i in ncol(far_covs)) {
  #    aux = abs(outer(far_covs[i, ], far_covs[i, ], FUN = "-"))
  #    temp = aux[lower.tri(aux)]
  #    cols_far_pairs = c(cols_far_pairs, which(temp<far_pairs[i]))
  #    rows_far_pairs = c(rows_far_pairs, row_count+(1:length(which(temp<far_pairs[i]))))
  #    vals_far_pairs = c(vals_far_pairs, rep(1, length(which(temp<far_pairs[i]))))
  #    row_count = max(rows_far_pairs)
  #  }  
  #}
  
  #! Near constraints
  rows_ind_near_pairs = list()
  if (!is.null(near_covs)) {
    rows_near = NULL
    cols_near = NULL
    vals_near = NULL
    n_near_covs = ncol(near_covs)
    for (j in 1:n_near_covs) {
      near_cov = near_covs[,j]
      #! Near on average constraints
      if (!is.null(near_groups)) {
        near_group = near_groups[j]
        row_ind_near_all = rep(row_count+1, n_dec)    
        col_ind_near_all = rep(1:n_dec, 1)			
        i_ind = rep(1:(n_tot-1), (n_tot-1):1)
        aux = matrix(rep(1:n_tot, n_tot), nrow = n_tot, byrow = F)
        j_ind = aux[lower.tri(aux)]
        vals_near_all = 	near_cov[i_ind]-near_cov[j_ind]-(near_group*rep(1, n_dec))	
        row_count	= max(row_ind_near_all)
      }
      #! Near on all pairs constraints
      if (!is.null(near_pairs)) {
        near_pair = near_pairs[j]
        aux = abs(outer(near_cov, near_cov, FUN = "-"))
        temp = as.vector(matrix(t(aux)[lower.tri(aux)], nrow = 1, byrow = TRUE))
        cols_ind_near_pairs = which(temp>near_pair)
        if (length(cols_ind_near_pairs)>0) {
          rows_ind_near_pairs[[j]] = row_count+(1:length(cols_ind_near_pairs))
          vals_near_pairs = rep(1, length(cols_ind_near_pairs))
          row_count	= max(rows_ind_near_pairs[[j]])
        }
        if (length(cols_ind_near_pairs)==0) {
          cols_ind_near_pairs = NULL
          rows_ind_near_pairs[[j]] = -1
          vals_near_pairs = NULL
        }
      }		
      #! Put together
      if (!is.null(near_groups) && is.null(near_pairs)) {
        rows_near = c(rows_near, row_ind_near_all)
        cols_near = c(cols_near, col_ind_near_all)
        vals_near = c(vals_near, vals_near_all)
      }
      if (is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        rows_near = c(rows_near, rows_ind_near_pairs[[j]])
        cols_near = c(cols_near, cols_ind_near_pairs)
        vals_near = c(vals_near, vals_near_pairs)
      }
      if (!is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] != -1) {
        rows_near = c(rows_near, row_ind_near_all, rows_ind_near_pairs[[j]])
        cols_near = c(cols_near, col_ind_near_all, cols_ind_near_pairs)
        vals_near = c(vals_near, vals_near_all, vals_near_pairs)
      }
      if (!is.null(near_groups) && !is.null(near_pairs) && rows_ind_near_pairs[[j]] == -1) {
        rows_near = c(rows_near, row_ind_near_all)
        cols_near = c(cols_near, col_ind_near_all)
        vals_near = c(vals_near, vals_near_all)
      }
    }
  }
  
  #! Near on average constraints
  #if (!is.null(near_covs)) {
  #  rows_near_all = rep(row_count+1, n_dec)
  #  cols_near_all = 1:n_dec
  #  i_ind = rep(1:(n_tot-1), (n_tot-1):1)
  #  aux = matrix(rep(1:n_tot, n_tot), nrow = n_tot, byrow = F)
  #  j_ind = aux[lower.tri(aux)]
  #  vals_near_all = near_covs[i_ind]-near_covs[j_ind]-(near_groups*rep(1, n_dec))
  #  row_count = max(rows_near_all)
  #}
  
  #! Near on all pairs constraints
  #if (!is.null(near_covs)) {
  #  aux = abs(outer(near_covs, near_covs, FUN = "-"))
  #  temp = aux[lower.tri(aux)]
  #  cols_near_pairs = which(temp>near_pairs)
  #  rows_near_pairs = row_count+(1:length(cols_near_pairs))
  #  vals_near_pairs = rep(1, length(cols_near_pairs))
  #  row_count = max(rows_near_pairs)
  #}
  #! Near on all pairs constraints
  #if (!is.null(near_covs)) {
  #  for (i in ncol(near_covs)) {
  #    aux = abs(outer(near_covs[i, ], near_covs[i, ], FUN = "-"))
  #    temp = aux[lower.tri(aux)]
  #    cols_near_pairs = c(cols_near_pairs, which(temp>near_pairs[i]))
  #    rows_near_pairs = c(rows_near_pairs, row_count+(1:length(which(temp<near_pairs[i]))))
  #    vals_near_pairs = c(vals_near_pairs, rep(1, length(which(temp<near_pairs[i]))))
  #    row_count = max(rows_near_pairs)
  #  }  
  #}
  
  #! Moment constraints
  if (!is.null(mom_covs)) {
    rows_mom_1 = NA
    cols_mom_1 = NA
    vals_mom_1 = NA
    rows_mom_2 = NA
    cols_mom_2 = NA
    vals_mom_2 = NA
    n_covs_m = ncol(mom_covs)
    for (i in 1:n_covs_m) {
      cov_m = mom_covs[, i]
      rows_mom_1 = c(rows_mom_1, rep(row_count+i, n_dec))
      cols_mom_1 = c(cols_mom_1, 1:n_dec)
      i_ind = rep(1:(n_tot-1), (n_tot-1):1)
      aux = matrix(rep(1:n_tot, n_tot), nrow = n_tot, byrow = F)
      j_ind = aux[lower.tri(aux)]
      vals_mom_1 = c(vals_mom_1, cov_m[i_ind]-cov_m[j_ind]-(mom_tols[i]*rep(1, n_dec)))
    }
    rows_mom_1 = rows_mom_1[-1]
    cols_mom_1 = cols_mom_1[-1]
    vals_mom_1 = vals_mom_1[-1]
    row_count = max(rows_mom_1)
    for (i in 1:n_covs_m) {
      cov_m = mom_covs[, i]
      rows_mom_2 = c(rows_mom_2, rep(row_count+i, n_dec))
      cols_mom_2 = c(cols_mom_2, 1:n_dec)
      i_ind = rep(1:(n_tot-1), (n_tot-1):1)
      aux = matrix(rep(1:n_tot, n_tot), nrow = n_tot, byrow = F)
      j_ind = aux[lower.tri(aux)]
      vals_mom_2 = c(vals_mom_2, cov_m[j_ind]-cov_m[i_ind]-(mom_tols[i]*rep(1, n_dec)))
    }
    rows_mom_2 = rows_mom_2[-1]
    cols_mom_2 = cols_mom_2[-1]
    vals_mom_2 = vals_mom_2[-1]
    rows_mom = c(rows_mom_1, rows_mom_2)
    cols_mom = c(cols_mom_1, cols_mom_2)
    vals_mom = c(vals_mom_1, vals_mom_2)
    row_count = max(rows_mom)
  }
  
  #! Exact matching constraints
  rows_exact = numeric()
  cols_exact = numeric()
  vals_exact = numeric()
  if (!is.null(exact_covs)) {
    n_exact_cats = ncol(exact_covs)
    for (i in 1:n_exact_cats) {
      rows_exact = c(rows_exact, rep(row_count+i, n_dec))
      cols_exact = c(cols_exact, 1:n_dec)
      dist_exact_cov = abs(outer(exact_covs[, i], exact_covs[, i], "-"))
      vals_exact = c(vals_exact, dist_exact_cov[lower.tri(dist_exact_cov)])
    }
    row_count	= max(rows_exact)
  }	
  
  #! Near-exact matching constraints
  rows_near_exact = numeric()
  cols_near_exact = numeric()
  vals_near_exact = numeric()
  if (!is.null(near_exact_covs)) {
    n_near_exact_cats = ncol(near_exact_covs)
    for (i in 1:n_near_exact_cats) {
      rows_near_exact = c(rows_near_exact, rep(row_count+j, n_dec))
      cols_near_exact = c(cols_near_exact, 1:n_dec)
      dist_near_exact_cov = abs(outer(near_exact_covs[, i], near_exact_covs[, i], "-"))
      vals_near_exact = c(vals_near_exact, dist_near_exact_cov[lower.tri(dist_near_exact_cov)])
    }
    row_count	= max(rows_near_exact)
  }
  
  #! Fine balance constraints
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
    j = 1
    for (i in 1:n_fine_cats) {
      rows_fine = c(rows_fine, rep(row_count+j, n_dec))
      cols_fine = c(cols_fine, 1:n_dec)
      dist_fine_cov = outer(fine_covs_2[, i], fine_covs_2[, i], "-")
      dist_fine_cov = t(dist_fine_cov)
      vals_fine = c(vals_fine, dist_fine_cov[lower.tri(dist_fine_cov)])
      if (j == 1) {
        rows_fine = rows_fine[-1]
        cols_fine = cols_fine[-1]
        vals_fine = vals_fine[-1]
      }
      j = j+1
    }
    row_count = max(rows_fine)
  }
  
  #! Near fine balance constraints
  if (!is.null(near_fine_covs)) {
    #! Transform near_fine_covs to a matrix of binary inds. for each cat. of each fine balancing covariate
    near_fine_covs_2 = rep(NA, nrow(near_fine_covs))
    n_near_fine_covs = ncol(near_fine_covs)
    j = 1
    for (i in 1:n_near_fine_covs) {
      aux = factor(near_fine_covs[, i])
      near_fine_covs_2 = cbind(near_fine_covs_2, diag(nlevels(aux))[aux,])
      if (j == 1) {
        near_fine_covs_2 = near_fine_covs_2[, -1]
      }
      j = j+1
    }
    n_near_fine_cats = ncol(near_fine_covs_2)
    j = 1
    for (i in 1:n_near_fine_cats) {
      for(h in 1:2) {
        rows_near_fine = c(rows_near_fine, rep(row_count+j, n_dec))
        cols_near_fine = c(cols_near_fine, 1:n_dec)
        dist_near_fine_cov = outer(near_fine_covs_2[, i], near_fine_covs_2[, i], "-")
        dist_near_fine_cov = t(dist_near_fine_cov)
        vals_near_fine = c(vals_near_fine, dist_near_fine_cov[lower.tri(dist_near_fine_cov)])
        if (j == 1) {
          rows_near_fine = rows_near_fine[-1]
          cols_near_fine = cols_near_fine[-1]
          vals_near_fine = vals_near_fine[-1]
        }
        j = j+1
      }
    }
    row_count = max(rows_near_fine)
  }
  
  #! n constraints
  if (!is.null(total_pairs)) {
    rows_n = rep(row_count+1, n_dec)
    cols_n = 1:n_dec
    vals_n = rep(1, n_dec)
    row_count = max(rows_n)
  }				
  
  #! Put together
  rows = c(rows_nbm, rows_far, rows_near, 
           rows_mom, rows_exact, rows_near_exact, rows_fine, rows_near_fine, rows_n)
  cols = c(cols_nbm, cols_far, cols_near, 
           cols_mom, cols_exact, cols_near_exact, cols_fine, cols_near_fine, cols_n)
  vals = c(vals_nbm, vals_far, vals_near, 
           vals_mom, vals_exact, vals_near_exact, vals_fine, vals_near_fine, vals_n)
  aux = cbind(rows, cols, vals)[order(cols), ]
  cnstrn_mat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  Amat = cnstrn_mat
  
  #! bvec
  bvec = rep(1,  length(table(rows_nbm)))
  #bvec = c(bvec, rep(0, length(table(rows_far_all))))
  #bvec = c(bvec, rep(0, length(table(rows_far_pairs))))
  
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
  
  #bvec = c(bvec, rep(0, length(table(rows_near_all))))
  #bvec = c(bvec, rep(0, length(table(rows_near_pairs))))
  
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
  
  bvec = c(bvec, rep(0, length(table(rows_mom))))
  if (!is.null(exact_covs)) {
    bvec = c(bvec, rep(0, ncol(exact_covs)))
  }
  if (!is.null(near_exact_covs)) {
    bvec = c(bvec, near_exact_devs)
  }
  bvec = c(bvec, rep(0, length(table(rows_fine))))
  #! near-fine
  if (!is.null(near_fine_covs)) {
    bvec_8_aux = rep(NA, length(rows_near_fine))
    bvec_8_aux[seq(1, length(rows_near_fine), 2)] = -near_fine_devs
    bvec_8_aux[seq(2, length(rows_near_fine), 2)] = near_fine_devs
    bvec = c(bvec, bvec_8_aux)
  }
  # total pairs
  if (!is.null(total_pairs)) {
    bvec = c(bvec, total_pairs)
  }
  
  
  #! ub
  ub = rep(1, n_dec)
  
  # sense
  sense = rep("L", length(table(rows_nbm)))
  #sense = c(sense, rep("G", length(table(rows_far_all))))
  #sense = c(sense, rep("E", length(table(rows_far_pairs))))
  
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
  
  #sense = c(sense, rep("L", length(table(rows_near_all))))
  #sense = c(sense, rep("E", length(table(rows_near_pairs))))
  
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
  
  sense = c(sense, rep("L", length(table(rows_mom))))
  if (!is.null(exact_covs)) {
    sense = c(sense, rep("E", ncol(exact_covs)))
  }
  if (!is.null(near_exact_covs)) {
    sense = c(sense, rep("L", ncol(near_exact_covs)))
  }
  sense = c(sense, rep("E", length(table(rows_fine))))
  sense = c(sense, rep(c("G", "L"), length(table(rows_near_fine))/2))
  sense = c(sense, rep("E", length(total_pairs)))
  
  # var_type
  if (approximate == 1) {
    var_type = rep("C", n_dec)
  }
  else {
    var_type = rep("B", n_dec)
  }	
  
  # Solve
  if (solver=="cplex") {
    #library("Rcplex")
    if (requireNamespace('Rcplex', quietly = TRUE)) {
      cat(format("  CPLEX optimizer is open..."), "\n")
      ptm = proc.time()
      out = Rcplex::Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = var_type,
                           control = list(trace = trace, round = round_cplex, tilim = t_max))
      time = (proc.time()-ptm)[3]
      
      # Output
      if (out$status==108) {
        cat(format("  Error: time limit exceeded, no integer solution!"), "\n")
        obj_val = NA
        obj_dist_mat = NA
        id_1 = NA
        id_2 = NA
        group_id = NA
        time = NA
      } else if (is.na(out$obj)) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj_val = NA
        obj_dist_mat = NA
        id_1 = NA
        id_2 = NA
        group_id = NA
        time = NA
      }
      
      if (!is.na(out$obj)) {
        cat(format("  Optimal matches found"), "\n")
        
        if (approximate == 1) {
          rel = .relaxation_n(n_tot, out$xopt, dist_mat, subset_weight, "cplex")
          out$xopt = rel$sol
          out$obj = rel$obj
          time = time + rel$time
        }
        
        i_ind = rep(1:(n_tot-1), (n_tot-1):1)
        aux = matrix(1:n_tot, nrow = n_tot, ncol = n_tot)
        j_ind = aux[lower.tri(aux)]
        
        group_1 = i_ind[out$xopt==1]
        group_2 = j_ind[out$xopt==1]
        max_groups = apply(cbind(group_1, group_2), 1, max)
        
        id_1 = group_1[max_groups<=n_tot]
        id_2 = group_2[max_groups<=n_tot]
        
        #! Group identifier
        group_id_1 = 1:(length(id_1))
        group_id_2 = 1:(length(id_2))
        group_id = c(group_id_1, group_id_2)
        
        obj_val = out$obj
        obj_dist_mat = sum(t(dist_mat)[lower.tri(dist_mat)] * out$xopt)
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver=="gurobi") {
    #library("gurobi")
    if (requireNamespace('gurobi', quietly = TRUE)) {
      cat(format("  Gurobi optimizer is open..."), "\n")
      model = list()
      model$modelsense = 'min'
      model$obj = cvec
      model$A = Amat
      model$sense = rep(NA, length(sense))
      model$sense[sense=="E"] = '='
      model$sense[sense=="L"] = '<='
      model$sense[sense=="G"] = '>='
      model$rhs = bvec
      model$vtypes = var_type
      model$ub = ub
      
      t_lim = list(TimeLimit = t_max, OutputFlag = trace)
      
      ptm = proc.time()
      out = gurobi::gurobi(model, t_lim)
      time = (proc.time()-ptm)[3]
      
      # Output
      if (out$status == "INFEASIBLE") {
        cat(format("  Error: problem infeasible!"), "\n")
        obj_val = NA
        obj_dist_mat = NA
        id_1 = NA
        id_2 = NA
        group_id = NA
        time = NA
      }
      
      if (out$status ==  "OPTIMAL" || out$status == "TIME_LIMIT") {
        if (out$status == "OPTIMAL") {
          cat(format("  Optimal matches found"), "\n")
        }
        else {
          cat(format("  Time limit reached, best suboptimal solution given"), "\n")
        }
        
        if (approximate == 1) {
          rel = .relaxation_n(n_tot, out$x, dist_mat, subset_weight, "gurobi")
          out$x = rel$sol
          out$objval = rel$obj
          time = time + rel$time
        }
        
        i_ind = rep(1:(n_tot-1), (n_tot-1):1)  
        aux = matrix(1:n_tot, nrow = n_tot, ncol = n_tot)
        j_ind = aux[lower.tri(aux)]
        
        group_1 = i_ind[out$x==1]
        group_2 = j_ind[out$x==1]
        max_groups = apply(cbind(group_1, group_2), 1, max)
        
        id_1 = group_1[max_groups<=n_tot]
        id_2 = group_2[max_groups<=n_tot]
        
        #! Group identifier
        group_id_1 = 1:(length(id_1))
        group_id_2 = 1:(length(id_2))
        group_id = c(group_id_1, group_id_2)
        
        obj_val = out$objval
        obj_dist_mat = sum(t(dist_mat)[lower.tri(dist_mat)] * out$x)
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver=="glpk") {
    #library("Rglpk")
    cat(format("  GLPK optimizer is open..."), "\n")
    dir = rep(NA, length(sense))
    dir[sense=="E"] = '=='
    dir[sense=="L"] = '<='
    dir[sense=="G"] = '>='
    bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                 upper = list(ind=c(1:length(ub)), val=ub))
    ptm = proc.time()
    out = Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = var_type, max = FALSE)
    time = (proc.time()-ptm)[3]
    
    # Output
    if (out$status!=0) {
      cat(format("  Error: problem infeasible!"), "\n")
      obj_val = NA
      obj_dist_mat = NA
      id_1 = NA
      id_2 = NA
      group_id = NA
      time = NA
    }
    
    if (out$status==0) {
      cat(format("  Optimal matches found"), "\n")
      
      if (approximate == 1) {
        rel = .relaxation_n(n_tot, out$solution, dist_mat, subset_weight, "glpk")
        out$solution = rel$sol
        out$optimum = rel$obj
        time = time + rel$time
      }
      
      i_ind = rep(1:(n_tot-1), (n_tot-1):1)
      aux = matrix(1:n_tot, nrow = n_tot, ncol = n_tot)
      j_ind = aux[lower.tri(aux)]
      
      group_1 = i_ind[out$solution==1]
      group_2 = j_ind[out$solution==1]
      max_groups = apply(cbind(group_1, group_2), 1, max)
      
      id_1 = group_1[max_groups<=n_tot]
      id_2 = group_2[max_groups<=n_tot]
      
      #! Group identifier
      group_id_1 = 1:(length(id_1))
      group_id_2 = 1:(length(id_2))
      group_id = c(group_id_1, group_id_2)
      
      obj_val = out$optimum
      obj_dist_mat = sum(t(dist_mat)[lower.tri(dist_mat)] * out$solution)
    }
  }
  
  if (solver=="symphony") {
    #library("Rsymphony")
    if (requireNamespace('Rsymphony', quietly = TRUE)) {
      cat(format("  Symphony optimizer is open..."), "\n")
      dir = rep(NA, length(sense))
      dir[sense=="E"] = '=='
      dir[sense=="L"] = '<='
      dir[sense=="G"] = '>='
      bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                   upper = list(ind=c(1:length(ub)), val=ub))
      ptm = proc.time()
      out = Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = var_type, max = FALSE, time_limit = t_max)
      time = (proc.time()-ptm)[3]
      
      # Output
      if (out$status==228) {
        cat(format("  Error: problem exceeded the time limit and no feasible solution is found!"), "\n")
        obj_val = NA
        obj_dist_mat = NA
        id_1 = NA
        id_2 = NA
        group_id = NA
        time = NA
      }
      else if (out$status!=0) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj_val = NA
        obj_dist_mat = NA
        id_1 = NA
        id_2 = NA
        group_id = NA
        time = NA
      }
      
      if (out$status==0) {
        cat(format("  Optimal matches found"), "\n")
        
        if (approximate == 1) {
          rel = .relaxation_n(n_tot, out$solution, dist_mat, subset_weight, "symphony")
          out$solution = rel$sol
          out$objval = rel$obj
          time = time + rel$time
        }
        
        i_ind = rep(1:(n_tot-1), (n_tot-1):1)
        aux = matrix(1:n_tot, nrow = n_tot, ncol = n_tot)
        j_ind = aux[lower.tri(aux)]
        
        group_1 = i_ind[out$solution==1]
        group_2 = j_ind[out$solution==1]
        max_groups = apply(cbind(group_1, group_2), 1, max)
        
        id_1 = group_1[max_groups<=n_tot]
        id_2 = group_2[max_groups<=n_tot]
        
        #! Group identifier
        group_id_1 = 1:(length(id_1))
        group_id_2 = 1:(length(id_2))
        group_id = c(group_id_1, group_id_2)
        
        obj_val = out$objval
        obj_dist_mat = sum(t(dist_mat)[lower.tri(dist_mat)] * out$solution)
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  return = list(obj_total = obj_val, obj_dist_mat = obj_dist_mat, id_1 = id_1, id_2 = id_2, 
                group_id = group_id,time = time)
}
