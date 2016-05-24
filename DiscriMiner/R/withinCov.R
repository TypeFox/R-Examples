#' Within-class Covariance Matrix
#' 
#' Calculates the within-class covariance matrix
#' 
#' When \code{div_by_n=TRUE} the covariance matrices are divided by n (number
#' of observations), otherwise they are divided by n-1
#' 
#' @param variables matrix or data frame with explanatory variables (No missing
#' values are allowed)
#' @param group vector or factor with group memberships (No missing values are
#' allowed)
#' @param div_by_n logical indicating division by number of observations
#' @author Gaston Sanchez
#' @seealso \code{\link{withinSS}}, \code{\link{betweenCov}},
#' \code{\link{totalCov}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load iris dataset
#'   data(iris)
#'   
#'   # within-class covariance matrix (dividing by n-1)
#'   withinCov(iris[,1:4], iris[,5])
#' 
#'   # within-class covariance matrix (dividing by n)
#'   withinCov(iris[,1:4], iris[,5], div_by_n=TRUE)
#'   }
#' 
withinCov <-
function(variables, group, div_by_n=FALSE)
{
  # within-class pooled covariance matrix
  # variables: matrix or data frame with explanatory variables
  # group: vector or factor with group memberships
  # div_by_n: logical indicating division by num of observations
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  X = verify_Xy$X
  y = verify_Xy$y
  
  # how many observations
  n = nrow(X)
  # how many variables
  p = ncol(X)
  # group levels and number of levels
  glevs = levels(y)
  ng = nlevels(y)
  # within cov matrix
  Within = matrix(0, p, p)
  for (k in 1:ng)
  {
    tmp <- y == glevs[k]
    nk = sum(tmp)
    if (div_by_n) {
      Wk = ((nk-1)/n) * var(X[tmp,])
    } else {
      # R version / SPSS
      #Wk = ((nk-1)/(n-ng)) * var(X[tmp,])
      Wk = ((nk-1)/(n-1)) * var(X[tmp,])
    }
    Within = Within + Wk
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
