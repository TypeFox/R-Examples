###########################################################################
##                                                                       ##
## residuals.mat() - 'residuals' method for MAT models                   ##
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
residuals.mat <- function(object, k, weighted = FALSE, ...)
  {
    auto <- FALSE
    if(missing(k))
      {
        auto <- TRUE
        if(weighted)
          k <- which.min(object$weighted$rmse)
        else
          k <- which.min(object$standard$rmse)
      }
    if(weighted)
      res <- object$weighted$resid[k, ]
    else
      res <- object$standard$resid[k, ]
    retval <- list(residuals = res, k = k, weighted = weighted,
                   auto = auto)
    class(retval) <- "residuals.mat"
    return(retval)
  }

print.residuals.mat <- function(x,
                                digits = min(3, getOption("digits") - 3), ...)
  {
    k <- x$k
    cat("\n")
    writeLines(strwrap("Modern Analogue Technique Residuals", prefix = "\t"))
    cat("\n")
    cat(paste("No. of analogues (k) :", k, "\n"))
    cat(paste("User supplied k?     :", !x$auto, "\n"))
    cat(paste("Weighted analysis?   :", x$weighted, "\n\n"))
    print.default(x$residuals, digits = digits)
    invisible(x)
  }
