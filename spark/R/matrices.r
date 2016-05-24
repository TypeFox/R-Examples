
#' Spark lines for rows of a matrix.
#'
#' @inheritParams spark
#' @inheritParams spark.default
#' @param data The matrix to plot.
#' @param common_scale Whether to plot all rows on a common scale.
#'   If \code{FALSE}, then all rows will have their own scales.
#' @param ... Passed to \code{spark.default}. E.g. it can be
#'   used to mark the minimum and/or maximum values.
#' @return Character scalar containing the spark line.
#'
#' @method spark matrix
#' @family spark
#' @export
#' @examples
#' spark(volcano)
#'
#' spark(t(EuStockMarkets), width="auto", common_scale = FALSE)

spark.matrix <- function(data, width = c("data", "auto", "screen"),
                         print = TRUE, common_scale = TRUE, ...) {

  if (common_scale) data <- t(apply(data, 1, scale_y, range = range(data)))

  res <- apply(data, 1, spark, width = width, print = FALSE, ...)

  if (print) {
    cat(res, sep = "\n")
    invisible(res)
  } else {
    res
  }

}
