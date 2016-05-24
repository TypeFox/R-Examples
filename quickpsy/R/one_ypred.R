#' Calculate the predictions of the response variable for one condition
#' \code{one_ypred} calculate the predictions of the response variable for one
#' condition
#' @keywords internal
#' @export
one_ypred <- function(d,log, groups, averages, x, psyfunguesslapses) {
  if (length(groups) != 0) averages <- semi_join(averages, d, by = groups)
  xseq <- averages[[x]]
  yseq <- psyfunguesslapses(xseq, d$par)
  data.frame(x = xseq, ypred = yseq)
}



