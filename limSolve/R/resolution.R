
##==============================================================================
## resolution  : calculates the resolution of equations and variables
## Given the singular value decomposition s, or the input matrix
## calculates the resolution of the equations (rows) and of the variables (columns)
## of the matrix
##==============================================================================

resolution <- function (s, tol=sqrt(.Machine$double.eps)) {

  if (is.numeric(s)) s <- svd(s)
  solvable <- sum (s$d>tol*s$d[1])       # number of sufficiently large singular values

  ## The resolution of the equations
  resolutioneq <- diag(s$u[,1:solvable]%*%t(s$u[,1:solvable]))

  ## The resolution of the variables
  resolutionvar <- diag(s$v[,1:solvable]%*%t(s$v[,1:solvable]))

  return(list(row=resolutioneq,      # resolution of the rows  (equations)
              col=resolutionvar,     # resolution of the columns (variables)
              nsolvable=solvable))   # number of solvable unknowns - rank

}

