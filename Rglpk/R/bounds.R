## bounds of objective coefficients

## fixes the GLPK bound types given a data.frame with bounds
## GLP_FR 1 free variable
## GLP_LO 2 variable with lower bound
## GLP_UP 3 variable with upper bound
## GLP_DB 4 double-bounded variable
## GLP_FX 5 fixed variable
glp_fix_bound_type <-
function(x)
{
  if(!inherits(x,"bound_table"))
    stop("'x' is not of class 'bound_table'")
  x$type <- ifelse(is.finite(x$lower),
                   ifelse(is.finite(x$upper), 4L, 2L),
                   ifelse(is.finite(x$upper), 3L, 1L))
  x$type[x$upper == x$lower] <- 5L
  x
}

## TODO: should be a generic function providing methods for
## different representations (e.g., a matrix, list of vectors, ...)
##                   

## A generic function which allows to take different dense and sparse
## representations of bounds.

as.glp_bounds <- function(x, ...)
  UseMethod("as.glp_bounds")

## No default representation.
as.glp_bounds.default <- function(x)
 stop("There is no default method for bounds representations.")

## returns identity
as.glp_bounds.bound_table <- function(x, n)
  x

## list -> GLPK bounds representation
as.glp_bounds.list <- function(x, n)
  glp_bounds(x, n)
  
glp_bounds <- function(x, n)
{
  ## General input validation
  ##if(!is.list(x))
  ##  stop("Bounds have to be of type list")

  ## Initialize default matrix
  bound_table <-
      expand.grid(type = rep.int(2L, n), lower = 0.0, upper = Inf)
  class(bound_table) <- c("bound_table", class(bound_table))
  
  ## Lower bounds
  lower <- x$lower
  ## check for zero-length bounds
  if( !any(unlist(lapply(lower, length))) )
    lower <- NULL

  if(!is.null(lower)){
    ## input validation
    glp_bounds_check_sanity(lower, n)
    if(any(lower[[1L]] == Inf))
      stop("Lower bound cannot be 'Inf'")
    ## if everything is OK set new lower bounds
    bound_table[lower[[1L]], 2L] <- lower[[2L]]
  }

  ## Upper bounds
  upper <- x$upper
  ## check for zero-length bounds
  if( !any(unlist(lapply(upper, length))) )
    upper <- NULL
  
  if(!is.null(upper)){
    ## input validation
    glp_bounds_check_sanity(upper, n)
    if(any(upper[[1L]] == -Inf))
      stop("Upper bound cannot be '-Inf'")
    ## so far, the same as with lower bounds but in addition we have to be
    ## sure that upper bounds are greater than or equal to lower bounds
    if(any(bound_table[upper[[1L]], 2L] > upper[[2L]]))
      stop("Upper bounds have to be greater than or equal to lower bounds")
    bound_table[upper[[1L]], 3L] <- upper[[2L]]
  }

  ## Fix bound types
  out <- glp_fix_bound_type(bound_table)
  out
}
  
glp_bounds_check_sanity <-
function(x, n)
{
  if(!is.numeric(x[[1L]]))
    warning("Bound indices not numeric. Coercing to integers ...")
  x[[1L]] <- as.integer(x[[1L]])
  if(length(x[[1L]]) != length(x[[2L]]))
    stop("Length of bound indices must be equal to the length of the corresponding bound values.")
  if(any(duplicated(x[[1L]])))
    stop("Duplicated entries in bound indices found.")
  if((max(x[[1L]]) > n))
    stop("Bound indices must not exceed number of objective coefficients.")
}
