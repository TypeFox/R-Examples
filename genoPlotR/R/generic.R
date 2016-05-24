################################################################################
# S3 method: generic functions
################################################################################
## trim generic function
trim <- function(x, ...){
  UseMethod("trim")
}
# slightly weird: given a xlim, returns rows for which the first numeric column
# is > xlim[1] and the second numeric column is < xlim[2] 
trim.default <- function(x, xlim=NULL, ...){
  if (!is.numeric(xlim)) stop("xlim must be numeric")
  if (length(xlim) != 2) stop("xlim must be length 2")
  num_col <- which(sapply(x, is.numeric))
  if (length(num_col) < 2) stop ("x should have at least 2 numeric columns")
  x <- x[x[,num_col[1]] >= xlim[1] & x[,num_col[2]] <= xlim[2]]
  x
}

## reverse generic function
reverse <- function(x, ...){
  UseMethod("reverse")
}
# default method
reverse.default <- function(x, ...){
  num_col <- which(sapply(x, is.numeric))
  if (length(num_col) < 2) stop ("x should have at least 2 numeric columns")
  tmp <- -x[,num_col[2]]
  x[,num_col[2]] <- -x[,num_col[1]]
  x[,num_col[1]] <- tmp
  x
}
