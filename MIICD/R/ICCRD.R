#' ICCRD : interval censored competing risks data
#' @name ICCRD 
#' @description Interval censored competing risks data. A data frame with 150 observations. the columns are :
#' \itemize{
#' \item left -> lower bound of the interval
#' \item right -> upper bound of the interval
#' \item status -> cause of failure (1 or 2)
#' \item treatment -> treatment (1 or 2)
#' \item cov2 -> another covariate ( continuous )
#' }
#' @details This dataset is given for demonstration purpose. 2 causes of failure are given, only cause 1 is interval censored. Right censored observations are indicated by \code{0} in the \code{status} column. 
#' @examples head(ICCRD)
#' @docType data
#' @keywords dataset
NULL


