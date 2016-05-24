################
# Dataset documentation
#
#' @title Blood and other measurements in diabetics [Hastie and Efron (2012)]
#'
#' @description
#' A repackaged diabetes dataset [Hastie and Efron (2012)] is a list of two different design  matrices 
#' and a response vector with 442 observations [Efron et. al. (2004)]
#' 
#' The x matrix has been standardized to have unit L2 norm in each
#' column and zero mean. The matrix x2 consists of x plus 54 interaction terms.
#'
#' @references Efron B., Hastie T., Johnstone I., and Tibshirani R. 
#' "Least Angle Regression" \emph{Annals of Statistics} 32:407-499, 2004.
#' @references Hastie T. and Efron B. (2012). lars: Least Angle Regression, Lasso and Forward Stagewise. 
#' R package version 1.1. http://CRAN.R-project.org/package=lars

#' @docType data
#' @keywords datasets
#' @name diabetes
#' @usage diabetes
#' 
#' @format A list of 3 data objects, 
#' \itemize{
#' \item x: A data frame with 10 variables and 442 observations
#' \item y: a numeric response vector of 442 observations
#' \item x2: a design matrix including interaction terms with 64 columns and 442 observations.
#' }
################
NULL
