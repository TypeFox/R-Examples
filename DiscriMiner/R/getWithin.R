#' Within-class Covariance Matrix
#' 
#' Calculates the estimated within-class covariance matrix
#' 
#' The obtained matrix is the estimated within-class covariance matrix (i.e.
#' within-class covariance matrix divided by its degrees of freedom \code{n-k},
#' where \code{n} is the number of observations and \code{k} is the number of
#' groups)
#' 
#' @param variables matrix or data frame with explanatory variables (No missing
#' values are allowed)
#' @param group vector or factor with group memberships (No missing values are
#' allowed)
#' @author Gaston Sanchez
#' @seealso \code{\link{withinCov}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load iris dataset
#'   data(iris)
#'   
#'   # estimated within-class covariance matrix (dividing by n-k)
#'   getWithin(iris[,1:4], iris[,5])
#' 
#'   # compared to the within-class covariance matrix (dividing by n-1)
#'   withinCov(iris[,1:4], iris[,5])
#'   }
#' 
getWithin <- 
function(variables, group)
{
  # within-class pooled covariance matrix
  # variables: matrix or data frame with explanatory variables
  # group: vector or factor with group memberships
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  X = verify_Xy$X
  y = verify_Xy$y
  
  Within = my_withinCov(X, y)
  Within
}
