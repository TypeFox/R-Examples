## the R-ported GNU Linear Programming kit
## solve function --- C Interface

Rglpk_solve_LP <-
function(obj, mat, dir, rhs, bounds = NULL, types = NULL, max = FALSE,
          control = list(), ...)
{
  ## validate direction of optimization
  if(!identical( max, TRUE ) && !identical( max, FALSE ))
    stop("'Argument 'max' must be either TRUE or FALSE.")
  direction_of_optimization <- as.integer(max)

  ## validate control list
  dots <- list(...)
  control[names(dots)] <- dots
  control <- .check_control_parameters( control )
  verb <- control$verbose

  ## match direction of constraints
  n_of_constraints <- length(dir)
  ## match relational operators to requested input
  direction_of_constraints <- match( dir, c("<", "<=", ">", ">=", "==") )

  if( any(is.na(direction_of_constraints)) )
    stop("Argument 'dir' must be either '<', '<=', '>', '>=' or '=='.")

  ## we need to verify that obj is a numeric vector
  ## FIXME: always use STMs?
  if(slam::is.simple_triplet_matrix(obj))
      obj <- as.matrix(obj)
  obj <- as.numeric(obj)
  n_of_objective_vars <- length( obj )

  constraint_matrix <- as.simple_triplet_matrix(mat)

  ## types of objective coefficients
  ## Default: "C"
  if(is.null(types))
    types <- "C"
  ## check if valid types
  if(any(is.na(match(types, c("I", "B", "C"), nomatch = NA))))
    stop("'types' must be either 'B', 'C' or 'I'.")
  ## replicate types to fit number of columns
  types <- rep(types, length.out = n_of_objective_vars)
  ## need a TRUE/FALSE integer/binary representation
  integers <- types == "I"
  binaries <- types == "B"

  ## do we have a mixed integer linear program?
  is_integer <- any( binaries | integers )

  ## bounds of objective coefficients
  bounds <- as.glp_bounds( as.list( bounds ), n_of_objective_vars )

  ## Sanity check: mat/dir/rhs
  if( !all(c(dim(mat)[ 1 ], length(rhs)) == n_of_constraints) )
      stop( "Arguments 'mat', 'dir', and/or 'rhs' not conformable." )
  ## Sanity check: mat, obj
  if( dim(mat)[ 2 ] != n_of_objective_vars )
      stop( "Arguments 'mat' and 'obj' not conformable." )

  ## call the C interface - this actually runs the solver
  x <- glp_call_interface(obj, n_of_objective_vars, constraint_matrix$i,
                          constraint_matrix$j, constraint_matrix$v,
                          length(constraint_matrix$v),
                          rhs, direction_of_constraints, n_of_constraints,
                          is_integer,
                          integers, binaries,
                          direction_of_optimization, bounds[, 1L],
                          bounds[, 2L], bounds[, 3L], verb)

  solution <- x$lp_objective_vars_values
  ## are integer variables really integers? better round values
  solution[integers | binaries] <-
    round( solution[integers | binaries])
  ## match status of solution
  status <- as.integer(x$lp_status)
  if(control$canonicalize_status){
      ## 0 -> optimal solution (5 in GLPK) else 1
      status <- as.integer(status != 5L)
  }
  list(optimum = sum(solution * obj), solution = solution, status = status)
}

## this function calls the C interface
glp_call_interface <-
function(lp_objective_coefficients, lp_n_of_objective_vars,
         lp_constraint_matrix_i, lp_constraint_matrix_j, lp_constraint_matrix_v,
         lp_n_of_values_in_constraint_matrix, lp_right_hand_side,
         lp_direction_of_constraints, lp_n_of_constraints, lp_is_integer,
         lp_objective_var_is_integer, lp_objective_var_is_binary,
         lp_direction_of_optimization,
         lp_bounds_type, lp_bounds_lower, lp_bounds_upper,
         verbose)
{
  out <- .C("R_glp_solve",
            lp_direction_of_optimization= as.integer(lp_direction_of_optimization),
            lp_n_of_constraints         = as.integer(lp_n_of_constraints),
            lp_direction_of_constraints = as.integer(lp_direction_of_constraints),
            lp_right_hand_side          = as.double(lp_right_hand_side),
            lp_n_of_objective_vars      = as.integer(lp_n_of_objective_vars),
            lp_objective_coefficients   = as.double(lp_objective_coefficients),
            lp_objective_var_is_integer = as.integer(lp_objective_var_is_integer),
            lp_objective_var_is_binary  = as.integer(lp_objective_var_is_binary),
            lp_is_integer               = as.integer(lp_is_integer),
            lp_n_of_values_in_constraint_matrix = as.integer(lp_n_of_values_in_constraint_matrix),
            lp_constraint_matrix_i      = as.integer(lp_constraint_matrix_i),
            lp_constraint_matrix_j      = as.integer(lp_constraint_matrix_j),
            lp_constraint_matrix_values = as.double(lp_constraint_matrix_v),
            lp_bounds_type              = as.integer(lp_bounds_type),
            lp_bounds_lower             = as.double(lp_bounds_lower),
            ## lp_n_of_bounds_l            = as.integer(length(lp_lower_bounds_i)),
            lp_bounds_upper             = as.double(lp_bounds_upper),
            ## lp_n_of_bounds_u            = as.integer(length(lp_upper_bounds_i)),
            lp_optimum                  = double(1),
            lp_objective_vars_values    = double(lp_n_of_objective_vars),
            lp_verbosity                = as.integer(verbose),
            lp_status                   = integer(1),
            NAOK = TRUE, PACKAGE = "Rglpk")
  out
}

## Convenience function for solving MILP objects
## upcoming ROI package (or for solving problems read with filereader)
.ROI_glpk_solve <- function(x, control = list()){
  if(!inherits(x, "MILP"))
    stop("'x' must be of class 'MILP'")
  if(is.null(control$verbose))
    control$verbose <- FALSE
  Rglpk_solve_LP(x$objective, x$constraints[[1]],
                 x$constraints[[2]], x$constraints[[3]],
                 types = x$types, max = x$maximum, bounds = x$bounds,
                 verbose = control$verbose)
}

