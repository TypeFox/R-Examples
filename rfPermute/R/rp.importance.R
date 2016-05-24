#' @title Extract rfPermute Importance Scores and p-values.
#' @description Extract a matrix of the observed importance scores
#'   and p-values from the object produced by a call to \code{rfPermute}
#' 
#' @param x An object produced by a call to \code{rfPermute}.
#' @param scale For permutation based measures, should the measures be divided 
#'   their "standard errors"?
#' @param sort.by character vector giving the importance metric(s) or p-values 
#'   to sort by. If \code{NULL}, defaults to \code{"MeanDecreaseAccuracy"} for 
#'   classification models and \code{"\%IncMSE"} for regression models.
#' @param decreasing logical. Should the sort order be increasing or decreasing?
#' 
#' @details p-values can be given to the \code{sort.by} argument by adding 
#'   '.pval' to the column name of the desired column from the \code{importance} 
#'   element of the \code{rfPermute} object.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{rfPermute}}
#' 
#' @export
#' 
rp.importance <- function(x, scale = TRUE, sort.by = NULL, decreasing = TRUE) {  
  if(!inherits(x, "rfPermute")) stop("'x' is not of class 'rfPermute'")
  if((!is.character(sort.by) & !is.vector(sort.by)) & !is.null(sort.by)) {
    stop("'sort.by' is not a character vector")
  }
  
  imp <- randomForest::importance(x, scale = scale)

  pval <- x$pval[, , if(scale) "scaled" else "unscaled"]
  colnames(pval) <- paste(colnames(pval), ".pval", sep = "")
  pred <- rownames(imp)
  vals <- do.call(cbind, lapply(1:ncol(imp), function(i) {
    cbind(imp[pred, i, drop = FALSE], pval[pred, i, drop = FALSE])
  }))
  
  if(is.null(sort.by)) sort.by <- ifelse(x$type == "regression", 
                                         "%IncMSE", "MeanDecreaseAccuracy")
  not.found <- sort.by[!(sort.by %in% colnames(vals))]
  if(length(not.found) > 0) {
    not.found <- paste(not.found, collapse = ", ")
    stop(paste("sort.by: ", not.found, " not found", sep = ""))
  }
  
  order.list <- lapply(sort.by, function(i) vals[, i, drop = FALSE])
  order.list <- c(order.list, decreasing = decreasing)
  imp <- vals[do.call(order, order.list), , drop = FALSE]
  if(!is.null(x$null.dist)) class(imp) <- c("rp.importance", class(imp))
  return(imp)
}