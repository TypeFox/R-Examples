#' Total Covariance Matrix
#' 
#' Calculates total covariance matrix
#' 
#' When \code{div_by_n=TRUE} the covariance matrices are divided by n (number
#' of observations), otherwise they are divided by n-1
#' 
#' @param variables matrix or data frame with explanatory variables (No missing
#' values are allowed)
#' @param div_by_n logical indicating division by number of observations
#' @author Gaston Sanchez
#' @seealso \code{\link{totalSS}}, \code{\link{betweenCov}},
#' \code{\link{withinCov}}
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load iris dataset
#'   data(iris)
#'   
#'   # total covariance matrix (dividing by n-1)
#'   totalCov(iris[,1:4])
#' 
#'   # total covariance matrix (dividing by n)
#'   totalCov(iris[,1:4], div_by_n=TRUE)
#'   }
#' 
totalCov <-
function(variables, div_by_n=FALSE)
{
  # Total covariance matrix
  # variables: matrix or data frame with explanatory variables
  # div_by_n: logical indicating division by num of observations
  
  # X matrix or data.frame
  if (!is.matrix(variables) && !is.data.frame(variables))
    stop("\nSorry, 'variables' must be a matrix")
  if (is.null(dim(variables)))
    stop("'variables' is not a matrix")
  # enforce variables as matrix
  if (!is.matrix(variables)) variables = as.matrix(variables)
  # no missing values allowed
  if (any(!is.finite(variables)))
    stop("infinite, NA or NaN values in 'variables'")
  # only numeric values
  if (!is.numeric(variables))
    stop("\nSorry, 'variables' must contain only numeric values")
  
  n = nrow(variables)
  Total = var(variables)
  if (div_by_n) {
    Total = ((n-1)/n) * Total
  }
    
  # add names
  if (is.null(colnames(variables))) {
    var_names = paste("X", 1:ncol(variables), sep="")
    dimnames(Total) = list(var_names, var_names)
  } else {
    dimnames(Total) = list(colnames(variables), colnames(variables))
  }
  # result  
  Total
}
