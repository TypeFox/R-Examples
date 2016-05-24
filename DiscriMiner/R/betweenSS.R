#' Between-class Sum of Squares Matrix
#' 
#' Calculates between-class sum of squares and cross product matrix (a.k.a.
#' between-class scatter matrix)
#' 
#' 
#' @param variables matrix or data frame with explanatory variables (No missing
#' values are allowed)
#' @param group vector or factor with group membership (No missing values are
#' allowed)
#' @author Gaston Sanchez
#' @seealso \code{\link{betweenCov}}, \code{\link{withinSS}},
#' \code{\link{totalSS}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load iris dataset
#'   data(iris)
#'   
#'   # between-class scatter matrix
#'   betweenSS(iris[,1:4], iris[,5])
#'   }
#' 
betweenSS <-
function(variables, group)
{
  # Between-class sum of squares and cross products matrix
  # variables: matrix or data frame with explanatory variables
  # group: vector or factor with group memberships
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  X = verify_Xy$X
  y = verify_Xy$y
  
  # group levels and number of levels
  glevs = levels(y)
  ng = nlevels(y)
  # global mean
  mean_all = colMeans(X)
  # matrix to store results
  Between = matrix(0, ncol(X), ncol(X))
  # calculate between Sum of squares
  for (k in 1:ng)
  {
    tmp <- y == glevs[k]
    nk = sum(tmp)
    mean_k = colMeans(X[tmp,])
    dif_k = mean_k - mean_all
    between_k = nk * dif_k %*% t(dif_k)
    Between = Between + between_k
  }
  # add names
  if (is.null(colnames(variables))) {
    var_names = paste("X", 1:ncol(X), sep="")
    dimnames(Between) = list(var_names, var_names)
  } else {
    dimnames(Between) = list(colnames(variables), colnames(variables))
  }
  # result
  Between
}
