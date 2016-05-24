## Reads linear programs from MPS files
## uses GLPK's facilities for reading these files
## Interface to GLPK's MPS file reader

## Since 2012-01-11: Rglpk supports MathProg files and retrieves
## variables names and constraints names. Thanks to Michael Kapler!

## Input: Path to a file specifying a mathematical program (MP),
##        the model specification language


## Output: an object of class 'MP_data_from_file' describing the MP

## Mathematical Programming (MP) Data Object (consisting of 2 parts)
## o MILP class
##$objective        ... a 1 x n simple_triplet_matrix representing objective coeffs
##$constraints[[1]] ... a m x n simple_triplet_matrix specifying the constraints
##$constraints[[2]] ... contains the m direction of constraints
##$constraints[[3]] ... vector of m right hand side values
##$bounds           ... a list with elements $upper and $lower. Each of which contains indices and bounds of objective variables
##$types            ... a character vector specifying which objective variable is of type 'binary' (B), continuous (C), or 'integer' (I)
##$maximum          ... can be either 'min' or 'max'
## o ATTRIBUTES
##$n_objective_vars ... number of objective variables
##$n_integer_vars   ... number of variables which are of type 'integer'
##$n_binary_vars    ... number of variables which are of type 'binary'
##$n_constraints    ... number of constraints
##$n_nonzero        ... number of values in constraint matrix
##$problem_name     ... name of the problem
##$objective_name   ... name of the problem
##$constraint_names ... names of the constraints
##$objective_vars_names ... names of the objective vars
##$file_name        ... absolute path to original data file
##$file_type        ... file type (currently 'MPS-fixed', 'MPS-free', 'CPLEX LP', 'MathProg')



Rglpk_read_file <- function(file, type = c("MPS_fixed", "MPS_free", "CPLEX_LP", "MathProg"), ignore_first_row = FALSE, verbose = FALSE){
  if(!file.exists(file))
    stop(paste("There is no file called", file, "!"))
  ## which file type to read from
  type <- match.arg(type)
  type_db <- c("MPS_fixed" = 1L,
               "MPS_free"  = 2L,
               "CPLEX_LP"  = 3L,
               "MathProg"  = 4L
               )
  type <- type_db[type]
  obj <- list(file = tools::file_path_as_absolute(file),
              type = type)
  ## read files in a two step approach
  ## first, retrieve meta data: e.g., the number of objective variables, etc.
  ## we need this to allocate memory on the R level for the result vectors
  meta_data <- glp_get_meta_data_from_file(obj, verbose)
  ## second, read all remaining data
  milp_data <- glp_retrieve_MP_from_file(meta_data, ignore_first_row, verbose)
  ## merge everything together
  MP_data <- glp_merge_MP_data(meta_data, milp_data)
  ## Post processing
  MP_data$type <- names(type_db[type_db == MP_data$type])
  ## and a direction '<='
  dir_db <- c("FR" = 1L, ">=" = 2L, "<=" = 3L, "DB" = 4L, "==" = 5L)
  MP_data$direction_of_constraints <- names(dir_db[MP_data$direction_of_constraints])
  ## we remove constraints marked as unbounded (GLP_FR)
  ind_FR <- MP_data$direction_of_constraints == "FR"
  if(any(ind_FR)){
      MP_data$constraint_matrix <- MP_data$constraint_matrix[ !ind_FR, ]
      MP_data$right_hand_side <- MP_data$right_hand_side[ !ind_FR ]
      MP_data$left_hand_side <- MP_data$left_hand_side[ !ind_FR ]
      MP_data$direction_of_constraints <- MP_data$direction_of_constraints[ !ind_FR ]
      MP_data$constraint_names <- MP_data$constraint_names[ !ind_FR ]
      MP_data$n_constraints <- MP_data$n_constraints - sum( ind_FR )
  }
  ## we add an additional constraint for double bounded constraints(GLP_DB)
  ind_DB <- MP_data$direction_of_constraints == "DB"
  if(any(ind_DB)){
      n_double_bounded <- sum( ind_DB )
      MP_data$constraint_matrix <- rbind( MP_data$constraint_matrix, MP_data$constraint_matrix[ ind_DB, ] )
      ## upper bounds already in rhs vector
      MP_data$direction_of_constraints[ ind_DB ] <- "<="
      length(MP_data$direction_of_constraints) <- length(MP_data$direction_of_constraints) + sum(ind_DB)
      MP_data$direction_of_constraints[is.na(MP_data$direction_of_constraints)] <- ">="
      length(MP_data$right_hand_side) <- length(MP_data$right_hand_side) + sum(ind_DB)
      MP_data$right_hand_side[is.na(MP_data$right_hand_side)] <- MP_data$left_hand_side[ind_DB ]
      MP_data$constraint_names <- c(MP_data$constraint_names, MP_data$constraint_names[ ind_DB ])
      MP_data$n_constraints <- MP_data$n_constraints + sum( ind_DB )
  }

  ## default is to have only continuous variables
  ## if any is binary or integer set the value accordingly
  types <- rep("C", length.out = MP_data$n_objective_vars)
  if(any(MP_data$objective_var_is_integer)){
    types[MP_data$objective_var_is_integer] <- "I"
  }
  if(any(MP_data$objective_var_is_binary)){
    types[MP_data$objective_var_is_binary] <- "B"
  }
  ## recalculate number of nonzeroes
  MP_data$n_values_in_constraint_matrix <- length(MP_data$constraint_matrix$v)

  ## build object we want to return
  ## First add MILP to the object
  out <- MILP(objective = MP_data$objective_coefficients,
              constraints = list(MP_data$constraint_matrix,
                                 MP_data$direction_of_constraints,
                                 MP_data$right_hand_side),
              bounds = MP_data$bounds,
              types = types,
              maximum = MP_data$maximize
              )
  attr(out, "n_objective_vars")     <- MP_data$n_objective_vars
  attr(out, "n_integer_vars")       <- MP_data$n_integer_vars
  attr(out, "n_binary_vars")        <- MP_data$n_binary_vars
  attr(out, "n_constraints")        <- MP_data$n_constraints
  attr(out, "n_nonzeros")           <- MP_data$n_values_in_constraint_matrix
  attr(out, "problem_name")         <- MP_data$problem_name
  attr(out, "objective_name")       <- MP_data$objective_name
  attr(out, "objective_vars_names") <- MP_data$objective_vars_names
  attr(out, "constraint_names")     <- MP_data$constraint_names
  attr(out, "file_type")            <- MP_data$type
  attr(out, "file_name")            <- MP_data$file


  class(out) <- c("MP_data_from_file", class(out))
  out
}

## First parse file to get some meta data of the LP/MILP
## (number of constraints/objective variables, direction of optimization, ...)
glp_get_meta_data_from_file <- function(x, verbose){
  res <- .C("Rglpk_read_file",
            file                          = as.character(x$file),
            type                          = as.integer(x$type),
            direction_of_optimization     = integer(1L),
            n_constraints                 = integer(1L),
            n_objective_vars              = integer(1L),
            n_values_in_constraint_matrix = integer(1L),
            n_integer_vars                = integer(1L),
            n_binary_vars                 = integer(1L),
            problem_name                  = character(1L),
            objective_name                = character(1L),
            verbosity                     = as.integer(verbose),
            PACKAGE = "Rglpk")
  ## free memory by deleting C-level problem object
  .C("Rglpk_delete_prob", PACKAGE = "Rglpk")
  res
}

## Retrieve all missing elements of the LP/MILP
glp_retrieve_MP_from_file <- function(x, ignore_first_row, verbose = FALSE){
  res <- .C("Rglpk_retrieve_MP_from_file",
            file                     = as.character(x$file),
            type                     = as.integer(x$type),
            n_constraints            = x$n_constraints,
            n_objective_vars         = x$n_objective_vars,
            objective_coefficients   = double(x$n_objective_vars),
            constraint_matrix_i      = integer(x$n_values_in_constraint_matrix),
            constraint_matrix_j      = integer(x$n_values_in_constraint_matrix),
            constraint_matrix_values = double(x$n_values_in_constraint_matrix),
            direction_of_constraints = integer(x$n_constraints),
            right_hand_side          = double(x$n_constraints),
            left_hand_side           = double(x$n_constraints),
            objective_var_is_integer = integer(x$n_objective_vars),
            objective_var_is_binary  = integer(x$n_objective_vars),
            bounds_type              = integer(x$n_objective_vars),
            bounds_lower             = double(x$n_objective_vars),
            bounds_upper             = double(x$n_objective_vars),
            lp_ignore_first_row      = as.integer(ignore_first_row),
            verbosity                = as.integer(verbose),
            constraint_names         = rep(character(1L), x$n_constraints),
            objective_vars_names     = rep(character(1L), x$n_objective_vars),
            PACKAGE = "Rglpk")
  ## free memory by deleting C-level problem object
  .C("Rglpk_delete_prob", PACKAGE = "Rglpk")
  ## lp_is_integer               = as.integer(lp_is_integer),

  ## replace infinity values
  res$bounds_lower <- replace(res$bounds_lower, res$bounds_lower == -.Machine$double.xmax, -Inf)
  res$bounds_upper <- replace(res$bounds_upper, res$bounds_upper == .Machine$double.xmax, Inf)
  ## in MPS definition first row is sometimes problematic. E.g., in MIPLIB2003
  ## it has to be removed!
  if(ignore_first_row){
    res$n_constraints <- res$n_constraints - 1
    ## zeros values in the constraint matrix have to be removed, these
    ## are the values from the first row
    to_remove <- which(res$constraint_matrix_values == 0)
    res$constraint_matrix_i <- res$constraint_matrix_i[-to_remove] - 1
    res$constraint_matrix_j <- res$constraint_matrix_j[-to_remove]
    res$constraint_matrix_values <- res$constraint_matrix_values[-to_remove]
    res$right_hand_side <- res$right_hand_side[-1]
    res$left_hand_side <- res$left_hand_side[-1]
    #res$direction_of_constraints <- res$direction_of_constraints[-length(res$right_hand_side)]
  }
  res
}

glp_merge_MP_data <- function(x, y){
  out <- list(objective_coefficients        = y$objective_coefficients,
              constraint_matrix             = simple_triplet_matrix(
                                                   y$constraint_matrix_i,
                                                   y$constraint_matrix_j,
                                                   y$constraint_matrix_values,
                                                   y$n_constraints,
                                                   y$n_objective_vars),
              direction_of_constraints      = y$direction_of_constraints,
              right_hand_side               = y$right_hand_side,
              left_hand_side                = y$left_hand_side,
              objective_var_is_integer      = as.logical(y$objective_var_is_integer),
              objective_var_is_binary       = as.logical(y$objective_var_is_binary),
              ## minimization if GLP_MIN (1L) or max if GLP_MAX (2L)
              maximize                      = x$direction_of_optimization == 2L,
              bounds                        = list(lower = list(ind = 1L:x$n_objective_vars,
                                                                val = y$bounds_lower),
                                                   upper = list(ind = 1L:x$n_objective_vars,
                                                                val = y$bounds_upper)),
              n_objective_vars              = x$n_objective_vars,
              n_integer_vars                = x$n_integer_vars,
              n_binary_vars                 = x$n_binary_vars,
              ## here from y because it might have changed -> ignore_first_row_parameter
              n_constraints                 = y$n_constraints,
              n_values_in_constraint_matrix = x$n_values_in_constraint_matrix,
              ## problem_name                  = x$problem_name,
              file                          = x$file,
              type                          = x$type,
              problem_name                  = x$problem_name,
              objective_name                = x$objective_name,
              constraint_names         = y$constraint_names,
              objective_vars_names     = y$objective_vars_names
              )
  out
}
