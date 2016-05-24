#' @aliases get.linProb,linProb,character-method
#' @rdname get.linProb-method
setMethod(f="get.linProb", signature=c("linProb", "character"),
  definition=function(object, type) {
    if ( !type %in% c("constraints", "direction", "rhs", "objective", "types", "bounds") ) {
      stop("get.cutList:: argument 'type' is not valid!\n")
    }
    if ( type == "constraints" ) {
      return(g_constraints(object))
    }
    if ( type == "direction" ) {
      return(g_direction(object))
    }
    if ( type == "rhs" ) {
      return(g_rhs(object))
    }
    if ( type == "objective" ) {
      return(g_objective(object))
    }
    if ( type == "types" ) {
      return(g_types(object))
    }
    if ( type == "bounds" ) {
      return(g_bounds(object))
    }
  }
)

#' @aliases set.linProb,linProb,character,list-method
#' @rdname set.linProb-method
setMethod(f="set.linProb", signature=c("linProb", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("objective", "direction", "rhs", "types",
        "removeCompleteConstraint", "addCompleteConstraint",
        "bounds", "constraints") ) {
      stop("set.linProb:: check argument 'type'!\n")
    }
    if ( type == "objective" ) {
      s_objective(object) <- input
    }
    if ( type == "direction" ) {
      s_direction(object) <- input
    }
    if ( type == "rhs" ) {
      s_rhs(object) <- input
    }
    if ( type == "types" ) {
      s_types(object) <- input
    }
    if  ( type == "removeCompleteConstraint" ) {
      s_remove_complete_constraint(object) <- input
    }
    if ( type == "addCompleteConstraint" ) {
      s_add_complete_constraint(object) <- input
    }
    if ( type == "bounds" ) {
      s_bounds(object) <- input
    }
    if ( type == "constraints" ) {
      s_constraints(object) <- input
    }
    validObject(object)
    return(object)
  }
)

#' @aliases calc.linProb,linProb,character,list-method
#' @rdname calc.linProb-method
setMethod(f="calc.linProb", signature=c("linProb", "character", "list"),
  definition=function(object, type, input) {
    if ( !type %in% c("solveProblem", "fixVariables") ) {
      stop("calc.linProb:: check argument 'type'!\n")
    }

    if ( type == "solveProblem" ) {
      return(c_solve_problem(object, input))
    }
    if ( type == "fixVariables" ) {
      return(c_fix_variables(object, input))
    }
  }
)

setMethod(f="g_constraints", signature=c("linProb"), definition=function(object) {
  return(object@constraints)
})

setMethod(f="g_direction", signature=c("linProb"), definition=function(object) {
  return(object@direction)
})

setMethod(f="g_rhs", signature=c("linProb"), definition=function(object) {
  return(object@rhs)
})

setMethod(f="g_objective", signature=c("linProb"), definition=function(object) {
  return(object@objective)
})

setMethod(f="g_types", signature=c("linProb"), definition=function(object) {
  return(object@types)
})

setMethod(f="g_bounds", signature=c("linProb"), definition=function(object) {
  return(list(upper=object@boundsUpper, lower=object@boundsLower))
})

setReplaceMethod(f="s_objective", signature=c("linProb", "list"), definition=function(object, value) {
  object@objective <- value[[1]]
  return(object)
})

setReplaceMethod(f="s_direction", signature=c("linProb", "list"), definition=function(object, value) {
  object@direction <- value[[1]]
  return(object)
})

setReplaceMethod(f="s_rhs", signature=c("linProb", "list"), definition=function(object, value) {
  object@rhs <- value[[1]]
  return(object)
})

setReplaceMethod(f="s_types", signature=c("linProb", "list"), definition=function(object, value) {
  object@types <- value[[1]]
  return(object)
})

setReplaceMethod(f="s_bounds", signature=c("linProb", "list"), definition=function(object, value) {
  # FIXME: check bounds input (lower|upper,...)
  object@boundsLower <- value$lower
  object@boundsUpper <- value$upper
  return(object)
})

setReplaceMethod(f="s_constraints", signature=c("linProb", "list"), definition=function(object, value) {
  object@constraints <- value[[1]]
  return(object)
})

setReplaceMethod(f="s_remove_complete_constraint", signature=c("linProb", "list"), definition=function(object, value) {
  input <- value[[1]]
  if ( !all(input %in% 1:length(g_rhs(object))) ) {
    stop("s_remove_complete_constraint:: elements of argument 'input' must be >=1 and <=",length(g_rhs(object)),"!\n")
  }
  object@constraints <- c_remove_row(g_constraints(object), input=list(input))
  object@direction <- g_direction(object)[-input]
  object@rhs <- g_rhs(object)[-input]
  return(object)
})

setReplaceMethod(f="s_add_complete_constraint", signature=c("linProb", "list"), definition=function(object, value) {
  input <- value[[1]]
  if ( g_nr_cols(g_constraints(object)) != g_nr_cols(g_constraints(input)) ) {
    stop("s_add_complete_constraint:: nrCols of 'object' and 'input' differ!\n")
  }
  if ( g_nr_constraints(input) > 0 ) {
    con <- g_constraints(input)
    for ( k in 1:g_nr_rows(con) ) {
      x <- g_row(con, input=list(k))
      object@constraints <- c_add_row(g_constraints(object), input=list(index=g_col_ind(x), values=g_values(x)))
    }
    object@direction <- c(g_direction(object), g_direction(input))
    object@rhs <- c(g_rhs(object), g_rhs(input))
  }
  return(object)
})

setMethod(f="c_solve_problem", signature=c("linProb", "list"), definition=function(object, input) {
  solver <- input[[1]]

  if ( !solver %in% c("glpk", "symphony", "lpSolve") ) {
    stop("'solver' needs to be eiter 'glpk', 'lpSolve' or 'symphony!\n'")
  }
  if ( solver == "glpk" ) {
    sol <- my.Rglpk_solve_LP(
      g_objective(object),
      g_constraints(object),
      g_direction(object),
      g_rhs(object),
      g_types(object),
      max = FALSE,
      bounds=g_bounds(object),
      verbose = FALSE)
  }
  if ( solver == "lpSolve" ) {
    stop("solving with 'lpSolve' not yet available!\n")
  }
  if ( solver == "symphony" ) {
    #sol <- Rsymphony_solve_LP(
    # g_objective(object),
    # g_constraints(object),
    # g_direction(object),
    # g_rhs(object),
    # bounds=g_bounds(object), #bounds
    # g_types(object),
    # max = FALSE)
    stop("solving with 'symphony' not yet available!\n")
  }
  if ( solver == "cplex" ) {
    #directionOrig <- g_direction(object),
    #sense <- rep(NA, length(directionOrig))
    #sense[directionOrig=="=="] <- "E"
    #sense[directionOrig=="<="] <- "L"
    #sense[directionOrig==">="] <- "G"
    #sol <- Rcplex(
    # g_objective(object),
    # g_constraints(object),
    # g_rhs(object),
    # Qmat = NULL,
    # lb = 0,
    # ub = 1,
    # control = list(),
    # objsense = "min",
    # sense = sense,
    # vtype = g_types(object),
    # n = 1)
    stop("solving with 'cplex' not yet available!\n")
  }
  return(sol)
})

setMethod(f="c_fix_variables", signature=c("linProb", "list"), definition=function(object, input) {
  lb <- input[[1]]
  ub <- input[[2]]
  primSupps <- input[[3]]

  if ( length(lb) != 1 | length(ub) != 1 ) {
    stop("c_fix_variables:: length of arguments 'lb' and 'ub' must equal 1!\n")
  }
  if ( !ub > lb ) {
    stop("c_fix_variables:: arguments 'ub' must be >= argument 'lb'!\n")
  }

  con <- g_constraints(object)
  rhs <- g_rhs(object)
  dir <- g_direction(object)
  obj <- g_objective(object)
  bounds <- g_bounds(object)
  nrVars <- g_nr_cols(con)

  my.lp <- make.lp(0, nrVars)

  set.objfn(my.lp, obj)
  set.bounds(my.lp, upper = bounds$upper$val)
  set.bounds(my.lp, lower = bounds$lower$val)

  for ( i in 1:g_nr_rows(con) ) {
    r <- g_row(con, input=list(i))
    cols <- g_col_ind(r)
    vals <- g_values(r)
    dd <- ifelse(dir[i]=="==", "=", dir[i])
    add.constraint(my.lp, vals, dd, rhs[i], indices=cols)
  }
  solve(my.lp)
  get.objective(my.lp)

  dual <- get.dual.solution(my.lp)
  dual <- dual[(2+length(get.rhs(my.lp))):length(dual)]
  if ( length(dual) != nrVars ) {
    stop("c_fix_variables:: length of arguments does not match!\n")
  }
  freqs <- obj
  sol <- get.variables(my.lp)
  sol[is.zero(sol)] <- 0
  sol[is.one(sol)] <- 1

  ### calculate reduced costs from dual solution
  reducedCosts <- freqs
  dualInd <- which(dual!=0)
  if ( length(dualInd) > 0 ) {
    reducedCosts[dualInd] <- sapply(dualInd, function(x) {
      min(freqs[x], dual[x])
    })
  }

  bas <- which(sol != 0)
  reducedCosts[bas] <- 0
  ### end calculation of reduced costs

  ### which variables could be set to zero?
  indSetToZero <- which(lb + reducedCosts >= ub)

  # do not fix primary suppressions to zero!
  indSetToZero <- indSetToZero[-which(indSetToZero %in% primSupps)]
  return(indSetToZero)
})

