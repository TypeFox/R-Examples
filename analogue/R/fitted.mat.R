###########################################################################
##                                                                       ##
## fitted.mat() - 'fitted' method for MAT models                         ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object on which method dispatch applied (Only 'mat')  ##
## k             - number of analogues to use. If missing 'k' is chosen  ##
##                 automatically as the 'k' that achieves lowest RMSE.   ##
## weighted      - Logical. Should the analysis use weighted mean of env ##
##                 data of analogues as fitted/estimated values?         ##
##                                                                       ##
###########################################################################
fitted.mat <- function(object, k, weighted = FALSE, ...)
  {
    auto <- FALSE
    if(missing(k))
      {
        auto <- TRUE
        if(weighted)
          k <- which.min(object$weighted$rmsep)
        else
          k <- which.min(object$standard$rmsep)
      }
    if(weighted) {
      ## need the row from this matrix
      est <- object$weighted$est[k,]
    } else {
      est <- object$standard$est[k,]
    }
    retval <- list(estimated = est, k = k, weighted = weighted,
                   auto = auto)
    class(retval) <- "fitted.mat"
    return(retval)
  }

print.fitted.mat <- function(x, digits = max(3, getOption("digits") - 3),
                             ...)
  {
    k <- x$k
    cat("\n")
    writeLines(strwrap("Modern Analogue Technique: Fitted values", prefix = "\t"))
    cat("\n")
    cat(paste("No. of analogues (k) :", k, "\n"))
    cat(paste("User supplied k?     :", !x$auto, "\n"))
    cat(paste("Weighted analysis?   :", x$weighted, "\n\n"))
    print.default(x$estimated, digits = digits)
    invisible(x)
  }
