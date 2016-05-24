## $Id$

print.ci2d <- function(x, ...)
  {
    cat("\n")
    cat("----------------------------\n")    
    cat("2-D Confidence Region Object\n")
    cat("----------------------------\n")    
    cat("\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
    cat("Number of data points: ", x$nobs, "\n")
    cat("Number of grid points: ", length(x$x), "x", length(x$y), "\n")
    cat("Number of confidence regions:", length(x$contours), "\n")
    cat("\n")

    tab <- data.frame(
                      "Region"=1:length(x$contours),
                      "CI Level"=as.numeric(names(x$contours)),
                      "X Min"=sapply(x$contours, function(XX) min(XX$x)),
                      "X Max"=sapply(x$contours, function(XX) max(XX$x)),
                      "Y Min"=sapply(x$contours, function(XX) min(XX$y)),
                      "Y Max"=sapply(x$contours, function(XX) max(XX$y))
                      )

    print(tab, row.names=FALSE, ...)

    x$summary <- tab

    class(x) <- c("ci2d.summary", "ci2d")

    return(x)
  }
