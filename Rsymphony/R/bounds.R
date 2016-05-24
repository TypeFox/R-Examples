## bounds of objective coefficients
## Rglpk bounding types for the time being ...

glp_fix_bound_type <- function(x){
  if(!inherits(x,"bound_table"))
    stop("'x' is not of class 'bound_table'")
  x$type <- ifelse(is.finite(x$lower),
                   ifelse(is.finite(x$upper), 4L, 3L),
                   ifelse(is.finite(x$upper), 2L, 1L))
  x$type[x$upper==x$lower] <- 5L
  x
}

## TODO: should be a generic function providing methods for
## different representations (e.g., a matrix, list of vectors, ...)
##                   
glp_bounds <- function(x, n){
  ## General input validation
  ##if(!is.list(x))
  ##  stop("Bounds have to be of type list")

  ## Initialize default matrix
  bound_table <- expand.grid(type=rep.int(2L,n), upper=0.0, lower=Inf)
  class(bound_table) <- c("bound_table", class(bound_table))
  
  ## Lower bounds
  lower <- x$lower
  if(!is.null(lower)){
    ## input validation
    glp_bounds_check_sanity(lower, n)
    if(any(lower[[1]]==Inf))
      stop("Lower bound cannot be 'Inf'")
    ## if everything is OK set new lower bounds
    bound_table[lower[[1]],2] <- lower[[2]]
  }

  ## Upper bounds
  upper <- x$upper
  if(!is.null(upper)){
    ## input validation
    glp_bounds_check_sanity(upper, n)
    if(any(upper[[1]]==-Inf))
      stop("Upper bound cannot be '-Inf'")
    ## so far, the same as with lower bounds but in addition we have to be
    ## sure that upper bounds are greater than or equal to lower bounds
    if(any(bound_table[upper[[1]],2] > upper[[2]]))
      stop("Upper bounds have to be greater than or equal to lower bounds")
    bound_table[upper[[1]],3] <- upper[[2]]
  }

  ## Fix bound types
  out <- glp_fix_bound_type(bound_table)
  out
}
  
glp_bounds_check_sanity <- function(x, n){
  if(!is.numeric(x[[1L]]))
    warning("Bound indices not numeric. Coercing to integers ...")
  x[[1L]] <- as.integer(x[[1L]])
  if(length(x[[1]]) != length(x[[2]]))
    stop("Length of bound indices must be equal to the length of the corresponding bound values!")
  if(any(duplicated(x[[1]])))
    stop("Duplicated entries in bound indices found!")
  if((max(x[[1]]) > n))
    stop("Bound indices must not exceed number of objective coefficients!")
}
