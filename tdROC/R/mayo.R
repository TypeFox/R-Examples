#' example data: Mayo Marker Data
#'
#' Two marker values with event time and censoring status for the subjects in Mayo PBC data
#' @docType data
#' @usage data(mayo)
#' @keywords datasets
#' @format A data frame with 312 observations and 4 variables: time (event time/censoring time), censor (censoring
#'        indicator), mayoscore4, mayoscore5. The two scores are derived from 4 and 5 covariates
#'        respectively.
#' @references Heagerty, P. J., & Zheng, Y. (2005). Survival model predictive accuracy and ROC curves. Biometrics, 61(1), 92-105.
#' @seealso \code{\link[survivalROC]{mayo}}
"mayo"
