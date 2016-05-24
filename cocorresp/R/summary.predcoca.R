"summary.predcoca" <- function(object, axes = c(1:min(6, object$n.axes)),
         display = c("species", "sites"), ...)
  {
    cocaScores <- scores(object, choices = axes, display = display)
    retval <- list(cocaScores = cocaScores, call = object$call,
                   lambda = object$lambda, namY = object$nam.dat$namY,
                   namX = object$nam.dat$namX,
                   loadings = lapply(object$loadings, `[`, , axes),
                   varianceExp = object$varianceExp,
                   totalVar = object$totalVar)
    class(retval) <- "summary.predcoca"
    retval
  }

