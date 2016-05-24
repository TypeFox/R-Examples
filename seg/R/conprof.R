# ------------------------------------------------------------------------------
# Function 'conprof'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
conprof <- function(data, grpID = 1, n = 999, graph = TRUE, add = FALSE, ...) {

  if (ncol(data) < 2 || !is.numeric(data))
    stop("'data' must be a numeric matrix with at least two columns", 
      call. = FALSE)
  if (length(grpID) > 1) {
    warning("'grpID' has more than one value", call. = FALSE)
    grpID <- grpID[1]
  }
  # if (missing(n))
  #   n <- formals(conprof.calc)$n
  val <- conprof.calc(data, grpID, n)
  
  if (graph) {
    if (!add) {
      plot(NA, xlim=c(0, 1.05), ylim=c(0, 1.05), xaxt="n", yaxt="n", ...)
      intrval <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
      axistxt <- intrval * 100
      axis(side = 1, at = intrval, labels = axistxt, ...)
      axis(side = 2, at = intrval, labels = axistxt, ...)      
    }
    
    lines(val$x, val$y, ...)
  }
  
  # Proportion of the group 
  p <- sum(data[,grpID]) / sum(data)
  
  above <- which(val$x >= p)
  below <- which(val$x < p)
  partA <- sum(val$y[above]) / n
  partB <- p - (sum(val$y[below]) / n)
  
  d <- (partB + partA) / (1 - p)
  list(x = val$x, y = val$y, d = d)
  
  # belowcurve <- sum(val$y) / n
  # abovecurve <- 1 - belowcurve
  # belowdotted <- (sum(x[,grpID]) / sum(x))
  # abovedotted <- 1 - belowdotted
}
