#' Correlation Ratio
#' 
#' Calculates the correlation ratio between a quantitaive variable and a
#' qualitative variable
#' 
#' No missing values are allowed
#' 
#' @param variable a single quantitative variable
#' @param group vector or factor with group memberships (qualitative variable)
#' @author Gaston Sanchez
#' @seealso \code{\link{FRatio}}, \code{\link{discPower}}
#' @references Tenenhaus, M. (2007) \emph{Statistique}. Dunod, Paris.
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # iris dataset 
#'   data(iris)
#'   
#'   # correlation ratio between Petal-Length and Species
#'   corRatio(iris$Petal.Length, iris$Species)
#'   }
#' 
corRatio <-
function(variable, group)
{
  # Correlation ratio
  # variable: vector with explanatory variable
  # group: vector or factor with group memberships
  
  if (!is.numeric(variable)) 
    stop("\nSorry, 'variable' must be a numeric vector")
  if (!is.factor(group)) group = as.factor(group)
  if (nlevels(group) == 1)
    stop("\nSorry, 'group' has only one category")
  # correlation ratio (eta)
  lm_temp = lm(variable ~ group)
  eta = summary(lm_temp)$r.squared
  # result
  eta
}
