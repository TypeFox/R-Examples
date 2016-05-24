#' Discriminant Power
#' 
#' Measures Discriminant Power of explanatory variables
#' 
#' No missing values are allowed
#' 
#' @param variables matrix or data frame with explanatory variables
#' @param group vector or factor with group membership
#' @return A data frame containing the following columns
#' @return \item{correl_ratio}{Correlation Ratios}
#' @return \item{wilks_lambda}{Wilks Lambda}
#' @return \item{F_statistic}{F-statistic}
#' @return \item{p_value}{p-value of F-statistic}
#' @author Gaston Sanchez
#' @seealso \code{\link{corRatio}}, \code{\link{FRatio}}
#' @references Tenenhaus M. (2007) \emph{Statistique}. Dunod, Paris.
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # bordeaux wines dataset
#'   data(bordeaux)
#' 
#'   # discriminant power
#'   dp = discPower(bordeaux[,2:5], bordeaux$quality)
#'   dp
#'   }
#' 
discPower <-
function(variables, group)
{
  # measure discriminant power of variables
  # variables: matrix or data frame with explanatory variables
  # group: vector or factor with group memberships
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  X = verify_Xy$X
  y = verify_Xy$y
  
  # how many observations
  n = nrow(X)
  # how many groups
  ng = nlevels(y)
  # how many variables
  p = ncol(X)
  # group levels and number of levels
  glevs = levels(y)
  
  # between-class covariance matrix
  B = my_betweenCov(X, y)
  # within-class covariance matrix
  W = matrix(0, p, p)
  for (k in 1:ng)
  {
    tmp <- y == glevs[k]
    nk = sum(tmp)
    Wk = ((nk-1)/(n-1)) * var(X[tmp,])
    W = W + Wk
  }
  
  # total covariance
  #V = ((n-1)/n) * var(X)
  V = var(X)
  
  ## Discriminant importance of explanatory variables
  # F-statistics
  F_statistic = ((n-ng)/(ng-1)) * (diag(B) / diag(W))
  p_value = 1 - pf(F_statistic, ng-1, n-ng)
  # Wilk's lambdas
  wilks_lambda = diag(W / V)
  # correlation ratios
  correl_ratio = diag(B) / diag(V)
  # table
  disc_power = data.frame(correl_ratio, wilks_lambda, F_statistic, p_value)
  disc_power
}
