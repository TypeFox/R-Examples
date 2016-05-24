#' F-Statistic Ratio
#' 
#' Calcualtes the F-statistic between a quantitative variable and a qualitative
#' variable
#' 
#' 
#' @param variable a quantitative variable
#' @param group a vector or factor with group memberships (i.e. qualitative
#' variable)
#' @return F-statistic and its p-value
#' @author Gaston Sanchez
#' @seealso \code{\link{discPower}}, \code{\link{corRatio}}
#' @references Tenenhaus M. (2007) \emph{Statistique}. Dunod, Paris.
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load bordeaux wines dataset
#'   data(bordeaux)
#' 
#'   # F-statistic ratio between temperature and quality
#'   FRatio(bordeaux$temperature, bordeaux$quality)
#'   }
#' 
FRatio <-
function(variable, group)
{
  # F ratio (anova)
  # variable: vector with explanatory variable
  # group: vector or factor with group memberships
  
  if (!is.numeric(variable)) 
    stop("\nSorry, argument 'x' must be a numeric vector")
  if (!is.factor(group)) group = as.factor(group)
  if (nlevels(group) == 1)
    stop("\nSorry, 'group' has only one category")
  # correlation ratio (eta)
  Ftest_temp = oneway.test(variable ~ group, var.equal=TRUE)
  fstat = Ftest_temp$statistic
  pval = Ftest_temp$p.value
  # result
  c(fstat, p_value=pval)
}
