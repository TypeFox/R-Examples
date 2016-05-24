###############################################################
## This function is adapted/modified based on the plot.cv
#    function from the glmnet package:
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via
#    Coordinate Descent.
## Journal of Statistical Software, 33(1), 1-22.
## URL http://www.jstatsoft.org/v33/i01/.
###############################################################

plot.cv.ernet <- function(x, sign.lambda = 1, ...) {
    cvobj <- x
    xlab <- "log(Lambda)"
    if (sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
    plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm, 
      ylim = range(cvobj$cvupper, cvobj$cvlower), xlab = xlab, 
      ylab = cvobj$name, type = "n")
    new.args <- list(...)
    if (length(new.args)) 
      plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper, 
      cvobj$cvlower, width = 0.01, col = "darkgrey")
    points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
    axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz), 
        tick = FALSE, line = 0)
    abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
    invisible()
} 
