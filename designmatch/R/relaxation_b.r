# Solves relaxation problem
.relaxation_b = function(n_t, n_c, coef, dist_mat, subset_weight, solver) {
  
  n_tot = n_t*n_c
  #! Part 1
  row_ind_1 = sort(rep(1:n_t, n_c))
  col_ind_1 = 1:n_tot
  val_1 = rep(1, n_tot)
  sense_1 = rep("L", n_t)
  bvec_1 = rep(1, n_t)
  
  row_ind_cur  = max(row_ind_1)
  
  #! Part 2
  row_ind_2 = rep(1:n_c, n_t)+row_ind_cur
  col_ind_2 = 1:n_tot
  val_2 = rep(1, n_tot)
  sense_2 = rep("L", n_c)
  bvec_2 = rep(1, n_c)
  
  row_ind = c(row_ind_1, row_ind_2)
  col_ind = c(col_ind_1, col_ind_2)
  vals = c(val_1, val_2)
  sense = c(sense_1, sense_2)
  bvec = c(bvec_1, bvec_2)
  
  aux = cbind(row_ind, col_ind, vals)[order(col_ind), ]
  Amat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  ub = Inf
  vtype = rep("B", n_tot)
  cvec = coef
  
  
  #### SOLVER PART #######
  
  if (solver == "cplex"){
    #library(Rcplex)
    if (requireNamespace('Rcplex', quietly = TRUE)) {
      ptm = proc.time()
      out = Rcplex::Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = vtype, n = 1,
                           control = list(trace = 0, round = 1), objsense = "max")
      time = (proc.time()-ptm)[3]
      
      if (out$status==108) {
        cat(format("  Error: time limit exceeded, no integer solution!"), "\n")
        obj = 0
        sol = NULL
      } else if (is.na(out$obj)) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      } else {
        sol = out$xopt
        if(is.null(dist_mat)) {
          obj = sum(-1*(rep(1, n_tot)) * out$xopt)
        } else {
          obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$xopt)
        }
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver == "gurobi") {
    #library(gurobi)
    if (requireNamespace('gurobi', quietly = TRUE)) {
      model = list()
      model$obj = cvec
      model$A = Amat
      model$sense = rep(NA, length(sense))
      model$sense[sense=="E"] = '='
      model$sense[sense=="L"] = '<='
      model$sense[sense=="G"] = '>='
      model$rhs = bvec
      model$vtypes = vtype
      model$ub = ub
      model$modelsense = "max"
      
      ptm = proc.time()
      out = gurobi::gurobi(model)
      time = (proc.time()-ptm)[3]
      
      if (out$status == "INFEASIBLE") {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      }
      
      if (out$status == "OPTIMAL") {
        sol = out$x
        if(is.null(dist_mat)) {
          obj = sum(-1*(rep(1, n_tot)) * out$x)
        } else {
          obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$x)
        }
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  if (solver == "symphony") {
    #library(Rsymphony)
    if (requireNamespace('Rsymphony', quietly = TRUE)) {
      dir = rep(NA, length(sense))
      dir[sense=="E"] = '=='
      dir[sense=="L"] = '<='
      dir[sense=="G"] = '>='
      
      bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                   upper = list(ind=c(1:length(ub)), val=ub))
      
      ptm = proc.time()
      out= Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
      time = (proc.time()-ptm)[3]
      
      if (out$status!=0) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj = 0
        sol = NULL
      } else {
        sol = out$solution
        if(is.null(dist_mat)) {
          obj = sum(-1*(rep(1, n_tot)) * out$solution)
        } else {
          obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$solution)
        }
      }
    } else {
      stop('suggested package not installed')
    }
  }
  
  # GLPK
  if (solver == "glpk") {
    #library(Rglpk)
    dir = rep(NA, length(sense))
    dir[sense=="E"] = '=='
    dir[sense=="L"] = '<='
    dir[sense=="G"] = '>='
    
    bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                 upper = list(ind=c(1:length(ub)), val=ub))
    
    ptm = proc.time()
    out= Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = TRUE)
    time = (proc.time()-ptm)[3]
    
    if (out$status!=0) {
      cat(format("  Error: problem infeasible!"), "\n")
      obj = 0
      sol = NULL
    } else {
      sol = out$solution
      if(is.null(dist_mat)) {
        obj = sum(-1*(rep(1, n_tot)) * out$solution)
      } else {
        obj = sum((as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))-(subset_weight*rep(1, n_t*n_c)))*out$solution)
      }
    }
    
  }
  
  return(list(sol = sol, obj = obj, time = time))
  
}