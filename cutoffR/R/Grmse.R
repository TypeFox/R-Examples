#' RMSE give imputed data matrix and the true matrix
#' @param ximp the imputed matrix
#' @param xtrue the true matrix
#' @return the RMSE
#' @export
#' @examples
#' data(hqmr.data)
Grmse <- function(ximp, xtrue) {
  err <- ximp - xtrue
  rmse <- sqrt((sum((unlist(err)^2), na.rm = TRUE))/(length(unlist(err)) - 
                                                          nmissing(err)))
  return(rmse)
}
