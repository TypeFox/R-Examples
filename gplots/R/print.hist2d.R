# $Id: hist2d.R 1471 2011-08-16 01:03:31Z warnes $

print.hist2d <- function(x, ...)
  {
    cat("\n")
    cat("----------------------------\n")    
    cat("2-D Histogram Object\n")
    cat("----------------------------\n")    
    cat("\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
    cat("Number of data points: ", x$nobs, "\n")
    cat("Number of grid bins: ", length(x$x), "x", length(x$y), "\n")
    cat("X range: (", min(x$x.breaks), ",", max(x$x.breaks), ")\n")
    cat("Y range: (", min(x$y.breaks), ",", max(x$y.breaks), ")\n")        
    cat("\n")
  }
