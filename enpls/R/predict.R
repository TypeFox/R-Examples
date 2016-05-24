#' Make Predictions from a Fitted Ensemble Partial Least Squares Model
#'
#' This function make predictions on new data by fitted enpls.en object.
#'
#' This function make predictions on new data by fitted enpls.en object.
#'
#' @param object An object of class \code{enpls.en}.
#' @param newx new data to predict with
#' @param ... Other parameters to be passed on to \code{predict}.
#'
#' @return A numeric vector containing the predicted values.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{enpls.en}} for performing ensemble PLS regression.
#'
#' @method predict enpls.en
#'
#' @export
#'
#' @references
#' Dongsheng Cao, Yizeng Liang, Qingsong Xu, Yifeng Yun, and Hongdong Li. 
#' "Toward better QSAR/QSPR modeling: simultaneous outlier detection and 
#' variable selection using distribution of model features." 
#' \emph{Journal of computer-aided molecular design} 25, no. 1 (2011): 67--80.
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' enpls.fit = enpls.en(x, y, MCtimes = 100)
#' y.pred = predict(enpls.fit, newx = x)
#' plot(y, y.pred, xlim = range(y), ylim = range(y))
#' abline(a = 0L, b = 1L)

predict.enpls.en = function(object, newx, ...) {

  if (missing(newx)) stop('Must provide newx')

  if (!inherits(object, 'enpls.en'))
    stop('This function only works for objects of class "enpls.en"')

  nmodel = length(object)

  predmat = matrix(NA, ncol = nmodel, nrow = nrow(newx))
  for (i in 1:nmodel) {
    predmat[, i] = predict(object[[i]][[1]], newx, object[[i]][[2]])
  }

  pred = rowMeans(predmat)

  return(pred)

}
