#' @title Visualizes the prevalence on the RBP curve.
#'
#' @description The prevalence is the proportion of a population having a specific condition.
#' In binary classification, the condition refers to whether the target variable has the value
#' \code{1}, that is, whether the target variable corresponds to the positive class.
#'
#' @template arg_obj
#' @template arg_plotvalues
#' @template arg_digits
#' @template arg_col
#' @template ret_invnull
#' @export
addPrevalence = function(obj, plot.values = TRUE, digits = 3L, col = "grey") {
  assertClass(obj, "RBPObj")

  # Compute 1-prevalence
  omp = obj$one.min.prev

  # Plot vertical lines where the distance between the lines reflects the prevalence
  abline(v = c(omp, 1), col = col)
  shape::Arrows(x0 = omp, x1 = 1L, y0 = -1L, y1 = -1L, 
    code = 3L, arr.adj = 1L, arr.col = col, col = col, lcol = col)

  # Should the value of the prevalence be plotted into the current plot?
  if (plot.values) {
    text(1 - (mean(obj$y) / 2), -1L, col = col,
      bquote(paste(hat(theta), " = ", .(round(obj$prev, digits)))), pos = 3L)
  } else {
    text(1 - (mean(obj$y) / 2), -1L, col = col, 
      expression(hat(theta)), pos = 3L)
  }

  return(invisible(NULL))
}

