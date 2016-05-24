
#                       Nodeplot: anti-boxplot
#                         Gilles Pujol 2006
# Modified by Frank Weber (2016): support model functions
# returning a matrix or a 3-dimensional array.

nodeplot <- function(x, xlim = NULL, ylim = NULL, labels = TRUE,
                     col = par("col"), pch = 21, bg = "white",
                     add = FALSE, at = NULL, y_col = NULL, y_dim3 = NULL, ...) {
  n <- nrow(x)
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1 : n
  }
  if (add) {
    par(new = TRUE)
  }
  
  # axes
  
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", type = "n", ...)
  if (class(labels) == "logical") {
    if (labels) {
      axis(side = 1, at = at, labels = rownames(x))
    } else {
      axis(side = 1, at = at, labels = FALSE, tick = FALSE)
    }
  } else if (class(labels) == "character") {
    axis(side = 1, at = at, labels = labels)
  }
  axis(side = 2)
  box()
  
  # bias
  
  if ("bias" %in% dimnames(x)[[2]]) {
    if(is.null(y_col) && is.null(y_dim3)){
      xx <- x[["original"]] - x[["bias"]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      xx <- x[, "original", y_col] - x[, "bias", y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      xx <- x[, "original", y_col, y_dim3] - x[, "bias", y_col, y_dim3]
    }
  } else {
    if(is.null(y_col) && is.null(y_dim3)){
      xx <- x[["original"]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      xx <- x[, y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      xx <- x[, y_col, y_dim3]
    }
  }
  
  # confidence intervals
  
  if (("min. c.i." %in% dimnames(x)[[2]]) & "max. c.i." %in% dimnames(x)[[2]]) {
    if(is.null(y_col) && is.null(y_dim3)){
      min_ci <- x[["min. c.i."]]
      max_ci <- x[["max. c.i."]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      min_ci <- x[, "min. c.i.", y_col]
      max_ci <- x[, "max. c.i.", y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      min_ci <- x[, "min. c.i.", y_col, y_dim3]
      max_ci <- x[, "max. c.i.", y_col, y_dim3]
    }
    for (i in 1 : n) {
      lines(c(at[i], at[i]), c(min_ci[i], max_ci[i]),
            col = col)
    }
  }
  
  # points
  
  points(at, xx, col = col, pch = pch, bg = bg)
}
