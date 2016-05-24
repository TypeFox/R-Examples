Rsymphony_solve_LP <-
function(obj, mat, dir, rhs, bounds = NULL, types = NULL, max = FALSE,
         verbosity = -2, time_limit = -1, node_limit = -1, gap_limit = -1,
         first_feasible = FALSE, write_lp = FALSE, write_mps = FALSE)
{
    ## Direction of optimization.
    if(!identical(max, TRUE) && !identical(max, FALSE))
        stop("'Argument 'max' must be either TRUE or FALSE.")
  
    nr <- nrow(mat)
    nc <- ncol(mat)

    ## Handle directions of constraints.
    TABLE <- c("L", "L", "E", "G", "G")
    names(TABLE) <- c('<', '<=', "==", ">", ">=")
    row_sense <- TABLE[dir]
    if(any(is.na(row_sense)))
        stop("Argument 'dir' must be one of '<', '<=', '>', '>=', or '=='.")
  
    ## Bounding support with using Rglpk bounds for the time being.
    bounds <- glp_bounds(as.list(bounds), nc)
    ## Use machine's max double values for infinities for the time being
    ## (as SYMPHONY does not know about IEEE 754 or C99 infinities).
    col_lb <-
        replace(bounds[, 2L], bounds[, 2L] == -Inf, -.Machine$double.xmax)
    col_ub <-
        replace(bounds[, 3L], bounds[, 3L] == Inf, .Machine$double.xmax)
  
    ## Note that the integer spec passed on is a vector of integer
    ## indicators, and that SYMPHONY has no native support for *binary*
    ## variables, so we pass treat these as integers <= 1.
    int <- if(is.null(types))
        logical(nc)
    else {
        if(!is.character(types) || !all(types %in% c("C", "I", "B")))
            stop("Invalid 'types' argument.")
        types <- rep(types, length.out = nc)
        col_ub[types == "B" & (col_ub > 1)] <- 1
        types != "C"
    }

    mat <- make_csc_matrix(mat)
  
    ## Call the C interface.
    out <- .C("R_symphony_solve",
              as.integer(nc),
              as.integer(nr),
              as.integer(mat$matbeg),
              as.integer(mat$matind),
              as.double(mat$values),
              as.double(col_lb),
              as.double(col_ub),
              as.integer(int),
              if(max) as.double(-obj) else as.double(obj),
              obj2 = double(nc),              
              as.character(paste(row_sense, collapse = "")),
              as.double(rhs),
              double(),
              objval = double(1L),
              solution = double(nc),
              status = integer(1L),
              verbosity = as.integer(verbosity),
              time_limit = as.integer(time_limit),
              node_limit = as.integer(node_limit),
              gap_limit = as.double(gap_limit),
              first_feasible = as.integer(first_feasible),
              write_lp = as.integer(write_lp),
              write_mps = as.integer(write_mps))

    ## Ensure that integer variables are really integer:
    solution <- out$solution
    solution[int] <- round(solution[int])
    
    status_db <- 
        c("TM_NO_PROBLEM" = 225L,
          "TM_NO_SOLUTION" = 226L,
          "TM_OPTIMAL_SOLUTION_FOUND" = 227L,
          "TM_TIME_LIMIT_EXCEEDED" = 228L,
          "TM_NODE_LIMIT_EXCEEDED" = 229L,
          "TM_ITERATION_LIMIT_EXCEEDED" = 230L,
          "TM_TARGET_GAP_ACHIEVED" = 231L,
          "TM_FOUND_FIRST_FEASIBLE" = 232L,
          "TM_FINISHED" = 233L,
          "TM_UNFINISHED" = 234L,
          "TM_FEASIBLE_SOLUTION_FOUND" = 235L,
          "TM_SIGNAL_CAUGHT" = 236L,
          "TM_UNBOUNDED" = 237L,
          "PREP_OPTIMAL_SOLUTION_FOUND" = 238L,
          "PREP_NO_SOLUTION" = 239L,
          "TM_ERROR__NO_BRANCHING_CANDIDATE" = -250L,
          "TM_ERROR__ILLEGAL_RETURN_CODE" = -251L,
          "TM_ERROR__NUMERICAL_INSTABILITY" = -252L,
          "TM_ERROR__COMM_ERROR" = -253L,
          "TM_ERROR__USER" = -275L,
          "PREP_ERROR" = -276L)
    status <- if(out$status == 227L)
        c("TM_OPTIMAL_SOLUTION_FOUND" = 0L)
    else
        status_db[match(out$status, status_db)]
    
    list(solution = solution,
         objval = sum(obj * solution),
         ## Equivalently,
         ##   if(max) - out$objval else out$objval
         status = status)
}        
