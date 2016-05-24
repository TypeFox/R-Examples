#! bmatch
bmatch = function(t_ind, dist_mat = NULL, subset_weight = NULL, n_controls = 1, total_pairs = NULL,
                   mom = NULL,
                   ks = NULL,
                   exact = NULL,
                   near_exact = NULL,
                   fine = NULL,
                   near_fine = NULL, 
                   near = NULL,
                   far = NULL,
                   #use_controls = NULL,
                   solver = NULL) {
  
  use_controls = NULL
  
  if (is.null(mom)) {
    mom_covs = NULL
    mom_tols = NULL
  } else {
    mom_covs = mom$covs
    mom_tols = mom$tols
  }
  
  if (is.null(ks)) {
    ks_covs = NULL
    ks_n_grid = 10
    ks_tols = NULL
  } else {
    ks_covs = ks$covs
    ks_n_grid = ks$n_grid
    ks_tols = ks$tols
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
  
  
  
  #! CALL ERROR HANDLING
  
  if (is.null(subset_weight)) {
    subset_weight = 0
  }
  
  #! Generate the parameters
  cat(format("  Building the matching problem..."), "\n")
  prmtrs = .problemparameters(t_ind, dist_mat, subset_weight, n_controls, total_pairs,
                             mom_covs, mom_tols,
                             ks_covs, ks_n_grid, ks_tols,
                             exact_covs,
                             near_exact_covs, near_exact_devs,
                             fine_covs,
                             near_fine_covs, near_fine_devs,
                             near_covs, near_pairs, near_groups,
                             far_covs, far_pairs, far_groups,
                             use_controls,
                             approximate)
  n_t = prmtrs$n_t
  n_c = prmtrs$n_c
  
  cvec = prmtrs$cvec
  Amat = prmtrs$Amat
  bvec = prmtrs$bvec
  ub = prmtrs$ub 
  sense = prmtrs$sense
  vtype = prmtrs$vtype
  c_index = prmtrs$c_index
  
  #! Find matches and calculate the elapsed time
  #! Gurobi
  if (solver == "gurobi") {
    #library(gurobi)
    if (requireNamespace('gurobi', quietly=TRUE)) {
      cat(format("  Gurobi optimizer is open..."), "\n")
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
      
      t_lim = list(TimeLimit = t_max, OutputFlag = trace)
      
      cat(format("  Finding the optimal matches..."), "\n")
      ptm = proc.time()
      out = gurobi::gurobi(model, t_lim)
      time = (proc.time()-ptm)[3]
      
      if (out$status == "INFEASIBLE") {
        cat(format("  Error: problem infeasible!"), "\n")
        obj_total = NA
        obj_dist_mat = NA
        t_id = NA
        c_id = NA
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
          rel = .relaxation_b(n_t, n_c, out$x, dist_mat, subset_weight, "gurobi")
          out$x = rel$sol
          out$objval = rel$obj
          time = time + rel$time
        }
        
        #! Matched units indexes
        t_id = sort(rep(1:n_t, n_c))[out$x==1]
        c_id = (c_index+n_t)[out$x==1]
        
        #! Group (or pair) identifier
        group_id_t = 1:(length(t_id))
        group_id_c = sort(rep(1:(length(t_id)), n_controls))
        group_id = c(group_id_t, group_id_c)
        
        #! Optimal value of the objective function
        obj_total = out$objval
        
        if (!is.null(dist_mat)) {
          obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$x))
        } else {
          obj_dist_mat = NULL
        }
      }
    } else {
      stop('suggested package not installed')
    }

  }
  
  #! CPLEX
  else if (solver == "cplex") {
    #library(Rcplex)
    if (requireNamespace('Rcplex', quietly=TRUE)) {
      cat(format("  CPLEX optimizer is open..."), "\n")
      cat(format("  Finding the optimal matches..."), "\n")
      ptm = proc.time()
      out = Rcplex::Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = vtype, n = 1, 
                           control = list(trace = trace, round = round_cplex, tilim = t_max))
      time = (proc.time()-ptm)[3]
      
      if (out$status==108) {
        cat(format("  Error: time limit exceeded, no integer solution!"), "\n")
        obj_total = NA
        obj_dist_mat = NA
        t_id = NA
        c_id = NA
        group_id = NA
        time = NA
      } else if (is.na(out$obj)) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj_total = NA
        obj_dist_mat = NA
        t_id = NA
        c_id = NA
        group_id = NA
        time = NA
      }
      
      if (!is.na(out$obj)) {
        cat(format("  Optimal matches found"), "\n")
        
        if (approximate == 1) {
          rel = .relaxation_b(n_t, n_c, out$xopt, dist_mat, subset_weight, "cplex")
          out$xopt = rel$sol
          out$obj = rel$obj
          time = time + rel$time
        }
        
        
        #! Matched controls indexes  
        t_id = sort(rep(1:n_t, n_c))[out$xopt==1]
        c_id = (c_index+n_t)[out$xopt==1]	
        
        #! Group (or pair) identifier
        group_id_t = 1:(length(t_id))
        group_id_c = sort(rep(1:(length(t_id)), n_controls))
        group_id = c(group_id_t, group_id_c)
        
        #! Optimal value of the objective function
        obj_total = out$obj
        
        if (!is.null(dist_mat)) {
          obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$xopt))
        } else {
          obj_dist_mat = NULL
        }
        
      }
    } else {
      stop('Suggested package not installed')
    }

    
  }
  
  #! GLPK
  else if (solver == "glpk") {
    #library(Rglpk)
    cat(format("  GLPK optimizer is open..."), "\n")
    dir = rep(NA, length(prmtrs$sense))
    dir[prmtrs$sense=="E"] = '=='
    dir[prmtrs$sense=="L"] = '<='
    dir[prmtrs$sense=="G"] = '>='
    
    bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                 upper = list(ind=c(1:length(ub)), val=ub))
    
    cat(format("  Finding the optimal matches..."), "\n")
    ptm = proc.time()
    out= Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = FALSE)
    time = (proc.time()-ptm)[3]
    
    if (out$status!=0) {
      cat(format("  Error: problem infeasible!"), "\n")
      obj_total = NA
      obj_dist_mat = NA
      t_id = NA
      c_id = NA
      group_id = NA
      time = NA
    }
    
    if (out$status==0) {
      cat(format("  Optimal matches found"), "\n")
      
      if (approximate == 1) {
        rel = .relaxation_b(n_t, n_c, out$solution, dist_mat, subset_weight, "glpk")
        out$solution = rel$sol
        out$optimum = rel$obj
        time = time + rel$time
      }
      
      
      #! Matched controls indexes	
      t_id = sort(rep(1:n_t, n_c))[out$solution==1]
      c_id = (c_index+n_t)[out$solution==1]	
      
      #! Group (or pair) identifier
      group_id_t = 1:(length(t_id))
      group_id_c = sort(rep(1:(length(t_id)), n_controls))
      group_id = c(group_id_t, group_id_c)
      
      #! Optimal value of the objective function
      obj_total = out$optimum
      
      if (!is.null(dist_mat)) {
        obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$solution))
      } else {
        obj_dist_mat = NULL
      }
      
    }
  }
  
  #! Symphony
  else {
    #library(Rsymphony)
    if (requireNamespace('Rsymphony', quietly=TRUE)) {
      cat(format("  Symphony optimizer is open..."), "\n")
      
      dir = rep(NA, length(prmtrs$sense))
      dir[prmtrs$sense=="E"] = '=='
      dir[prmtrs$sense=="L"] = '<='
      dir[prmtrs$sense=="G"] = '>='
      
      bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))),
                   upper = list(ind=c(1:length(ub)), val=ub))
      
      cat(format("  Finding the optimal matches..."), "\n")
      ptm = proc.time()
      out= Rsymphony::Rsymphony_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = vtype, max = FALSE, time_limit = t_max)
      time = (proc.time()-ptm)[3]
      
      if (out$status==228) {
        cat(format("  Error: problem exceeded the time limit and no feasible solution is found!"), "\n")
        obj_total = NA
        obj_dist_mat = NA
        t_id = NA
        c_id = NA
        group_id = NA
        time = NA
      }
      else if (out$status!=0) {
        cat(format("  Error: problem infeasible!"), "\n")
        obj_total = NA
        obj_dist_mat = NA
        t_id = NA
        c_id = NA
        group_id = NA
        time = NA
      }
      
      if (out$status==0) {
        cat(format("  Optimal matches found"), "\n")
        
        if (approximate == 1) {
          rel = .relaxation_b(n_t, n_c, out$solution, dist_mat, subset_weight, "symphony")
          out$solution = rel$sol
          out$objval = rel$obj
          time = time + rel$time
        }
        
        #! Matched controls indexes	
        t_id = sort(rep(1:n_t, n_c))[out$solution==1]
        c_id = (c_index+n_t)[out$solution==1]	
        
        #! Group (or pair) identifier
        group_id_t = 1:(length(t_id))
        group_id_c = sort(rep(1:(length(t_id)), n_controls))
        group_id = c(group_id_t, group_id_c)
        
        #! Optimal value of the objective function
        obj_total = out$objval
        
        if (!is.null(dist_mat)) {
          obj_dist_mat = sum(c(as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) * out$solution))
        } else {
          obj_dist_mat = NULL
        }
        
      }
    } else {
      stop('suggested package not installed')
    }

  }
  #! Output
  return(list(obj_total = obj_total, obj_dist_mat = obj_dist_mat, 
              t_id = t_id, c_id = c_id, group_id = group_id, time = time, status=out$status))
}