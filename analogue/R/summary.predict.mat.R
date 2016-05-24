###########################################################################
##                                                                       ##
## predict.mat() - 'predict' method for MAT models                       ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## object        - object on which method dispatch applied (Only 'mat')  ##
##                                                                       ##
###########################################################################
summary.predict.mat <- function(object, ...)
  {
    class(object) <- "summary.predict.mat"
    return(object)
  }

print.summary.predict.mat <- function(x,
                                      digits = max(3, getOption("digits") - 3),
                                      ...)
{
    print.predict.mat(x)
    if(!is.null(x$bootstrap)) {
        k.model <- x$model$k
        k.boot <- x$bootstrap$k
        dat <- data.frame(Obs = x$observed,
                          Est = x$model$estimated[k.model, ],
                          Resid = x$model$residuals[k.model, ],
                          Boot.Est = x$bootstrap$estimated[, k.boot],
                          Boot.Resid = x$bootstrap$residuals[, k.boot],
                          s1 = x$sample.errors$s1[, k.boot],
                          s2 = x$sample.errors$s2[, k.boot],
                          RMSEP = x$sample.errors$rmsep[, k.boot])
        cat("\nTraining set assessment:\n\n")
        print(dat, digits = digits)
    } else {
        k.model <- x$model$k
        dat <- data.frame(Obs = x$observed,
                          Est = as.numeric(x$model$estimated),
                          Resid = x$model$residuals$residuals)
        rownames(dat) <- names(x$model$estimated)
    }
    cat("\nTraining set assessment:\n\n")
    print(dat, digits = digits)
    invisible(x)
}
