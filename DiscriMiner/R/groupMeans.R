#' Group Means
#' 
#' Calculates means for each group
#' 
#' 
#' @param variables matrix or data frame with explanatory variables (may
#' contain missing values)
#' @param group vector or factor with group memberships
#' @param na.rm logical indicating whether missing values should be removed
#' @return matrix of group means (with variables in the rows, and groups in the
#' columns)
#' @author Gaston Sanchez
#' @seealso \code{\link{groupVars}}, \code{\link{groupStds}},
#' \code{\link{groupMedians}}, \code{\link{groupQuants}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # dataset iris
#'   data(iris)
#' 
#'   # group means
#'   groupMeans(iris[,1:4], iris[,5])
#'   }
#' 
groupMeans <-
function(variables, group, na.rm=FALSE)
{
  # Calculate means by group
  # variables: matrix or data frame with explanatory variables
  # group: vector or factor with group memberships
  # na.rm: logical indicating whether missing values should be removed
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=na.rm)
  X = verify_Xy$X
  y = verify_Xy$y

  # how many groups
  ng = nlevels(y)
  # matrix with group means
  Means = matrix(0, ncol(X), ng)
  for (j in 1:ncol(X))
  {
    Means[j,] = tapply(X[,j], y, FUN=mean, na.rm=na.rm)
  }
  # add names
  if (is.null(colnames(X))) {
    rownames(Means) = paste("X", 1:ncol(X), sep="")
  } else {
    rownames(Means) = colnames(X)
  }
  colnames(Means) = levels(y)
  # results
  Means
}
