#' Creates the curve for one condition
#' \code{one_curve} Creates the curve for one condition
#' @keywords internal
#' @export

one_curve <- function(d, xmin, xmax, log, groups, limits, psyfunguesslapses) {
  if (length(groups) != 0) limits <- semi_join(limits, d, by = groups)
  if (!is.null(xmin)) limits$xmin <- xmin
  if (!is.null(xmax)) limits$xmax <- xmax



  xseq <- seq(limits$xmin, limits$xmax, length = 300)
  yseq <- psyfunguesslapses(xseq, d$par)


  if (log) xseq <- exp(xseq)
  data.frame(x = xseq, y = yseq)
}
