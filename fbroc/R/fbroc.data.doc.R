#' Examples of predictions for ROC curve construction
#' 
#' Contains simulated data that can be used as examples for generating ROC curves. Both a continuous
#' and a discrete predictor are included. For both cases there is a version with outliers and one
#' without.
#' 
#' @format A data.frame with 160 rows and 5 variables:
#' \describe{
#'   \item{True.Class}{True class label of the observation}
#'   \item{Cont.Pred}{Predictions for which the binormal model for ROC curves holds. Predictions for
#'   both the positive and negative class follows a normal distribution with unit standard deviation
#'   and means 2 and 0 respectively.}
#'   \item{Cont.Pred.Outlier}{Same as above, with some extreme outliers in the negative class.}
#'   \item{Disc.Pred}{Example of a discrete predictor. Predictions for the negative class are integer 
#'   values between 1 and 8, positive samples have integer predictions between 7 and 14.}
#'   \item{Disc.Pred.Outlier}{Same as above, with some extreme outliers in the negative class.}
#' }
#' 
"roc.examples"