###########################################################################
##                                                                       ##
## summary.analog - summary method for class 'analog'                    ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object of class 'analog'                              ##
## display       - character vector listing which aspects of results to  ##
##                 display. One or more of 'dist', 'names' and           ##
##                 'quantiles', selecting the actual dissimilarities,    ##
##                 the names of analogues, and the quantiles of the      ##
##                 distribution of dissimilarities in modern training    ##
##                 set.                                                  ##
## k             - number of analogues to use.                           ##
## probs         - numeric vector, giving the probabilities of the       ##
##                 distribution to return quantiles for. See ?quantiles. ##
##                                                                       ##
###########################################################################
summary.analog <- function(object, display = c("dist", "names", "quantiles"),
                           k = 10, probs = c(0.01, 0.02, 0.05, 0.1, 0.2), ...)
  {
    summ <- list()
    display <- match.arg(display, several.ok = TRUE)
    #ord <- apply(object$analogs, 2, order)
    if("dist" %in% display)
      summ$dists <- apply(object$analogs, 2, function(x){ x[order(x)][1:k]})
    if("names" %in% display)
      {
        ord <- apply(object$analogs, 2, order)
        summ$names <- apply(ord, 2,
                            function(ord, x) {rownames(x)[ord][1:k]},
                            object$analogs)
      }
    if("quantiles" %in% display)
      summ$quantiles <- quantile(as.dist(object$train), probs = probs)
    summ$method <- object$method
    summ$call <- object$call #attr(object, "call")
    summ$k <- k
    class(summ) <- "summary.analog"
    summ
  }

print.summary.analog <- function(x, 
                                 digits = min(3, getOption("digits") - 4), ...)
  {
    method <- x$method
    .call <- deparse(x$call)
    closest <- x$k
    cat("\n")
    writeLines(strwrap("Analogue matching for fossil samples", prefix = "\t"))
    cat(paste("\nCall:", .call, "\n"))
    cat(paste("Dissimilarity:", method, "\n"))
    cat(paste("k-closest:", closest, "\n"))
    if(!is.null(x$quantiles))
      {
        cat("\nPercentiles of the dissimilarities for the training set:\n\n")
        print(x$quantiles, digits)
        cat("\n")
      }
    if(!is.null(x$names))
      {
        cat("k-closest analogues\n\n")
        tbl <- as.matrix(format(x$names, digits))
        tbl <- cbind(k = as.integer(1:closest), tbl)
        rownames(tbl) <- rep("", nrow(tbl))
        print(tbl, quote = FALSE, right = TRUE, print.gap = 2)
        cat("\n")
      }
    if(!is.null(x$dists))
      {
        cat("Dissimilarities for k-closest analogues\n\n")
        tbl <- as.matrix(format(x$dists, digits = digits))
        tbl <- cbind(k = as.integer(1:closest), tbl)
        rownames(tbl) <- rep("", nrow(tbl))
        print(tbl, quote = FALSE, right = TRUE, print.gap = 2)
        cat("\n")
      }
    invisible(x)
  }
