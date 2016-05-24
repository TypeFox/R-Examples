###########################################################################
##                                                                       ##
## bootstrap() - Generic function to calculate bootstrap statistics for  ##
##               transfer function models - only a method for class mat  ##
##               currently available                                     ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
###########################################################################
summary.bootstrap.mat <- function(object, ...)
  {
    class(object) <- "summary.bootstrap.mat"
    #class(object) <- "summary.bootstrap"
    return(object)
  }

print.summary.bootstrap.mat <- function(x,
#print.summary.bootstrap <- function(x,
                                        digits = max(3, getOption("digits") - 3),
                                        ...)
  {
    print.bootstrap.mat(x, digits = digits, ...)
    #print.bootstrap(x)
    cat("\nBootstrap estimated values for training set:\n")
    with(x$bootstrap, print(estimated[,k], digits = digits, ...))
    if(!is.null(x$predictions)) {
      cat(paste("\nPredicted values based on a",
                ifelse(x$weighted, " weighted", ""),
                " model with ", x$apparent$k,
                "-closest analogues\n", sep = ""))
      if(x$auto)
        cat("(k chosen from model with lowest RMSEP)\n\n")
      else
        cat("\n")
      with(x$predictions$model, print(predicted[k,], digits = digits, ...))
    }
    if(!is.null(x$predictions)) {
      cat(paste("\nBoostrap predicted values based on a",
                ifelse(x$weighted, " weighted", ""),
                " model with ", x$bootstrap$k,
                "-closest analogues\n", sep = ""))
      if(x$auto)
        cat("(k chosen from model with lowest RMSEP)\n\n")
      else
        cat("\n")
      with(x$predictions$bootstrap, print(predicted[,k], digits = digits, ...))
    }
    cat("\nTraining set assessment:\n\n")
    k.model <- x$model$k
    k.boot <- x$bootstrap$k
    dat <- data.frame(Obs = x$observed,
                      Est = x$model$estimated[k.model,],
                      Resid = x$model$residuals[k.model,],
                      Boot.Est = x$bootstrap$estimated[,k.boot],
                      Boot.Resid = x$bootstrap$residuals[,k.boot],
                      s1 = x$sample.errors$s1[,k.boot],
                      s2 = x$sample.errors$s2[,k.boot],
                      RMSEP = x$sample.errors$rmsep[,k.boot])
    print(dat, digits = digits, ...)
    invisible(x)
  }
