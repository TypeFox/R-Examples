#' @title Visualize the PEV on the RBP curve.
#'
#' @description
#' The PEV measure is the difference between the conditional expectation of the predicted
#' probabilities (conditional on the two groups that are determined by the target variable).
#' The PEV measure can be visually obtained by the RBP curve, namely by the difference of the
#' two areas that are Highlighted with \code{addPEV}.
#'
#' @template arg_obj
#' @template arg_plotvalues
#' @template arg_showinfo
#' @param text.col [\code{character(1)} | \code{numeric(1)}]\cr
#'   Text color, used when \code{plot.values = TRUE}, otherwise ignored.
#'   Default is \dQuote{black}.
#' @template arg_col
#' @template ret_invnull
#' @export
addPEV = function(obj, plot.values = TRUE, show.info = TRUE, 
  text.col = "black", col = rgb(0, 0, 0, 0.25)) {

  # Check arguments
  assertClass(obj, "RBPObj")
  assertFlag(plot.values)
  assert(checkString(text.col), checkNumeric(text.col))
  assertVector(col, len = 1L)

  # Store values of obj
  x1 = obj$axis.x
  y1 = obj$axis.y
  omp = obj$one.min.prev
  E0 = round(obj$e0, 4L)
  E1 = round(obj$e1, 4L)

  # Values of the x-axis that are greater than 1 - prevalence
  ind = x1 > omp
  
  # Highlight the area of the probabilities when Y=1
  polygon(x = c(x1[ind], 1, x1[ind][1L]), y = c(y1[ind], 1, 1),
    border = NA, col = col)

  # Highlight the area of the probabilities when Y=0
  polygon(x = c(x1[1L], x1[!ind], x1[length(x1[!ind])]),
    y = c(0, y1[!ind], 0), border = NA, col = col)

  # Add values for E1 and E0 into the plot
  if (plot.values) {
    text(min(x1), 0L, adj = 0:1, col = text.col,
      labels = bquote(paste(hat(E)[0], " = ", .(E0))))
    text(x1[length(x1[!ind])+1], 1L, adj = 0:1, col = text.col,
      labels = bquote(paste(hat(E)[1], " = ", .(E1))))
  }

  # Show message
  if (show.info) {
    messagef("E1: Mean predicted probabilities for Y=1: %s", E0)
    messagef("E0: Mean predicted probabilities for Y=0: %s", E1)
    messagef("PEV = E1 - E0: %s", round(obj$pev, 4L))
  }

  return(invisible(NULL))
}

