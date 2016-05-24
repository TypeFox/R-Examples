#' Total Sum of Squares Matrix
#' 
#' Calculates the total sum of squares and cross product matrix (a.k.a. total
#' scatter matrix)
#' 
#' 
#' @param variables matrix or data frame with explanatory variables
#' @author Gaston Sanchez
#' @seealso \code{\link{totalCov}}, \code{\link{betweenSS}},
#' \code{\link{withinSS}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load iris dataset
#'   data(iris)
#'   
#'   # total scatter matrix
#'   totalSS(iris[,1:4])
#'   }
#' 
totalSS <-
function(variables)
{
  # Total sum of squares matrix
  # variables: matrix or data frame with explanatory variables
  
  # X matrix or data.frame
  if (!is.matrix(variables) && !is.data.frame(variables))
    stop("\nSorry, 'variables' must be a matrix")
  if (is.null(dim(variables)))
    stop("'variables' is not a matrix")
  # enforce X as matrix
  if (!is.matrix(variables)) variables = as.matrix(variables)
  # no missing values allowed
  if (any(!is.finite(variables)))
    stop("infinite, NA or NaN values in 'variables'")
  # only numeric values
  if (!is.numeric(variables))
    stop("\nSorry, 'variables' must contain only numeric values")
  
  X = scale(variables, scale=FALSE)
  Total = t(X) %*% X
  
  # add names
  if (is.null(colnames(variables))) {
    var_names = paste("X", 1:ncol(X), sep="")
    dimnames(Total) = list(var_names, var_names)
  } else {
    dimnames(Total) = list(colnames(variables), colnames(variables))
  }
  # result
  Total
}
