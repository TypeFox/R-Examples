#' Within-class Sum of Squares Matrix
#' 
#' Calculates within-class sum of squares and cross product matrix (a.k.a.
#' within-class scatter matrix)
#' 
#' 
#' @param variables matrix or data frame with explanatory variables (No missing
#' values are allowed)
#' @param group vector or factor with group memberships (No missing values are
#' allowed)
#' @author Gaston Sanchez
#' @seealso \code{\link{withinCov}}, \code{\link{betweenSS}},
#' \code{\link{totalSS}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load iris dataset
#'   data(iris)
#'   
#'   # within-class scatter matrix
#'   withinSS(iris[,1:4], iris[,5])
#'   }
#' 
withinSS <-
function(variables, group)
{
  # Pooled within-class sum of squares and cross products
  # variables: matrix or data frame with explanatory variables
  # group: vector or factor with group memberships
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  X = verify_Xy$X
  y = verify_Xy$y
  
  # how many observations
  nrx = nrow(X)
  # how many variables
  ncx = ncol(X)
  # group levels and number of levels
  glevs = levels(y)
  ng = nlevels(y)
  # within cov matrix
  Within = matrix(0, ncx, ncx)
  for (k in 1:ng)
  {
    # select obs of k-th group
    tmp <- y == glevs[k]
    # mean k-th group
    mean_k = colMeans(X[tmp,])
    # center k-th group matrix
    Xk = scale(X[tmp,], center=mean_k, scale=FALSE)
    # add k-th intra-class SS matrix
    Within = Within + t(Xk) %*% Xk
  }
  # add names
  if (is.null(colnames(variables))) {
    var_names = paste("X", 1:ncol(X), sep="")
    dimnames(Within) = list(var_names, var_names)
  } else {
    dimnames(Within) = list(colnames(variables), colnames(variables))
  }
  # result
  Within
}
