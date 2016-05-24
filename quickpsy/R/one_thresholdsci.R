#' Calculates the confidence intervals for the threshold of one condition
#' \code{one_thresholdsci} calculates the confidence intervals for the threshold
#' of one condition
#' @keywords internal
#' @export
one_thresholdsci <- function(d, ci, method) {
  if (method == 'percent') {
    threinf <- quantile(d$thre, .5*(1 - ci))[[1]]
    thresup <- quantile(d$thre, 1 - .5*(1 - ci))[[1]]
  }
  data.frame(threinf, thresup)
}

