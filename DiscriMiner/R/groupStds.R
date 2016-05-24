#' Group Standard Deviations
#' 
#' Calculates the standard deviations for each group
#' 
#' 
#' @param variables matrix or data frame with explanatory variables (may
#' contain missing values)
#' @param group vector or factor with group memberships
#' @param na.rm logical indicating whether missing values should be removed
#' @return matrix of group standard deviations (with variables in the rows, and
#' groups in the columns)
#' @author Gaston Sanchez
#' @seealso \code{\link{groupMeans}}, \code{\link{groupVars}},
#' \code{\link{groupMedians}}, \code{\link{groupQuants}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # dataset iris
#'   data(iris)
#' 
#'   # group standard deviations
#'   groupStds(iris[,1:4], iris[,5])
#'   }
#' 
groupStds <-
function(variables, group, na.rm=FALSE)
{
  # Calculate std deviations by group
  # variables: matrix or data frame with explanatory variables
  # group: vector or factor with group memberships
  # na.rm: logical indicating whether missing values should be removed
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=na.rm)
  X = verify_Xy$X
  y = verify_Xy$y
  
  # how many groups
  ng = nlevels(y)
  # matrix with group std deviations
  Stds = matrix(0, ncol(X), ng)
  for (j in 1:ncol(X))
  {
    Stds[j,] = tapply(X[,j], y, FUN=sd, na.rm=na.rm)
  }
  # add names
  if (is.null(colnames(X))) {
    rownames(Stds) = paste("X", 1:ncol(X), sep="")
  } else {
    rownames(Stds) = colnames(X)
  }
  colnames(Stds) = levels(y)
  # results
  Stds
}
