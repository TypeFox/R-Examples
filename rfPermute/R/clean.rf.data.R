#' @title Clean Random Forest Input Data
#' @description Removes cases for a Random Forest classification model 
#'   with missing data and predictors that are constant.
#' 
#' @param x columns used as predictor variables as character or numeric vector.
#' @param y column used as response variable as character or numeric.
#' @param data data.frame containing \code{x} and \code{y} columns.
#' @param max.levels maximum number of levels in response variable \code{y}.
#' 
#' @return a data.frame containing cleaned data.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom stats complete.cases
#' @export
#' 
clean.rf.data <- function(x, y, data, max.levels = 30) {
  data <- as.data.frame(data)
  if (is.null(colnames(data))) colnames(data) <- 1:ncol(data)
  x <- setdiff(x, y)
  if (is.numeric(x)) x <- colnames(data)[x]
  if (is.numeric(y)) y <- colnames(data)[y]
  sub.df <- data[, c(y, x)]
  sub.df <- sub.df[complete.cases(sub.df), , drop = TRUE]
  
  delete.pred <- sapply(x, function(pred){
    pred.vec <- sub.df[[pred]]
    if (length(unique(pred.vec)) <= 1) return(pred)
    if (is.factor(pred.vec) & (nlevels(pred.vec) > max.levels)) return(pred)
    NULL
  })
  delete.pred <- unlist(delete.pred)
  delete.pred <- delete.pred[!is.null(delete.pred)]
  
  if (length(delete.pred) > 0) x <- setdiff(x, delete.pred)
  if (is.factor(sub.df[[y]]) & nlevels(sub.df[[y]][, drop = TRUE]) < 2) return(NULL)
  sub.df[, c(y, x)]
}
