#' @title Visualizes a measure for well calibration on the RBP curve.
#'
#' @description A measure for a well calibrated model can be obtained by
#' grouping the predicted probabilities via deciles yielding 10 groups.
#' The equally collored areas belong to a specific group. When each of the two  
#' equally collored areas are similar, the model is well calibrated.
#'
#' @template arg_obj
#' @template arg_plotvalues
#' @param subplot.control [\code{list}]
#'   A named list of arguments that will be passed to \code{\link{barplot}}.
#'   Additionally, you can set \code{diff = TRUE} to plot differences of the 
#'   equally collored areas or \code{diff = FALSE} to directly plot the areas 
#'   of the equally collored areas in juxtaposed bars.
#' @param col [\code{character} | \code{numeric}]\cr
#'   A specification for the the plotting color for the areas.
#' @param pos [\code{list}] 
#'   A named List that determines the \code{x} and \code{y} positioning of a subplot that
#'   compares the areas in additional barplots (see \code{\link{subplot}}).
#'   Can be \code{NA} for no additional subplot.
#'   Default is \code{pos = NULL} for an auto positioning in the topleft quadrant.
#' @return A matrix that contains the average of the \dQuote{probabilities within deciles}
#'   conditional on Y.
#' @export

addWellCalib = function(obj, plot.values = TRUE, subplot.control = list(diff = TRUE),
  col = shape::greycol(10L, interval = c(0.3, 1)), pos = NULL) {

  # Check arguments
  assertClass(obj, "RBPObj")
  assertFlag(plot.values)
  assertVector(col, min.len = 1L, max.len = 10L)
  if(is.null(pos)) pos = list(x =  c(0.15, 0.85 - obj$prev), y = c(0.3, 0.9))
  if(any(is.na(pos))) pos = list(x = as.numeric(NA), y = as.numeric(NA))
  assertList(pos, types = "numeric", len = 2)
  assertSubset(names(pos), c("x", "y"))
  
  if("diff"%in%names(subplot.control)){
    diff = subplot.control$diff
    subplot.control$diff = NULL
  } else diff = FALSE
  name = c("las", "xaxt", "col", "main", "height", "beside")
  if(any(name%in%names(subplot.control))){
    stopf("arguments %s not allowed in subplot.control", 
      paste(name, collapse=", "))
  }
  
  # Store values of obj
  pred = obj$pred
  y = obj$y
  n = obj$n
  x1 = obj$axis.x
  y1 = obj$axis.y
  eps = y - pred
  
  # Set range for deciles
  q = seq(0, 1, by = 0.1)
  up = 1 - q
  lo = -q

  # Matrix that will contain the average of the "probabilities within deciles" conditional on Y.
  areas = matrix(ncol = 2L, nrow = (length(q) - 1))
  row.names(areas) = sprintf("[%s, %s]", q[1:10], q[2:11])
  colnames(areas) = c("Y = 0", "Y = 1")
  
  for (i in 1:(length(q) - 1)) {
    # Residuals for different conditions
    eps1 = eps[pred < q[i + 1L] & pred >= q[i] & y == 1L]
    eps0 = eps[pred < q[i + 1L] & pred >= q[i] & y == 0L]
    
    areas[i, ] = c(sum(eps0), sum(eps1)) / n
    
    # Highlight the area for the probabilities bounded by the i-th and (i+1)-th decile
    if (length(eps0) != 0L) 
      areaCalib(x1, y1, lo[i], col[i], 
        label = round(areas[i, "Y = 0"], 4), plot.values = plot.values, pos = 4)
    if (length(eps1) != 0L) 
      areaCalib(x1, y1, up[i], col[i], 
        label = round(areas[i, "Y = 1"], 4), plot.values = plot.values, pos = 2)
  }
  
  # Add subplot
  if (!any(is.na(pos))) {    
    TeachingDemos::subplot(fun = {
      if(!diff){
        args = list(height = t(abs(areas)), beside = TRUE, col = rep(col, each = 2L), main = "Area")
      }else{
        args = list(height = rowSums(areas), col = col, main = "Area difference")
      }
      do.call("barplot", append(append(args, list(las = 2L, xaxt = "n")), subplot.control))
      }, x = pos$x, y = pos$y, pars = list(mar = c(0, 0, 1, 0) + 0.1)
    )
  }
  
  return(invisible(areas))
}

# helper for highlighting the area bounded by deciles of predicted risks
areaCalib = function(x1, y1, thres, col, label, plot.values, pos) {
  wg <- ifelse(pos == 4, 1, 3)
  ind = (y1 < thres) & (y1 >= thres - 0.1)
  if (sum(ind) != 0L) {
    polygon(x = c(min(x1[ind]), x1[ind], max(x1[ind])),
      y = c(0, y1[ind], 0), border = 1L, col = col)
    if (plot.values) {
      TeachingDemos::shadowtext(
        x = mean(c(rep(min(x1[ind]), pos-1), rep(max(x1[ind]), wg))),
        #(min(x1[ind])+max(x1[ind]))/2, 
        y = mean(c(thres, thres - 0.1)), bg = "white",
        labels = label, 
        pos = pos, col = col)
    }
  }
}
