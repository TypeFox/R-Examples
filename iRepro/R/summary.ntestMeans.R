summary.ntestMeans <-
function(object, ...){
  t <- as.data.frame(c(object$statistic,object$p.value))
  rownames(t) <- c("Chi-squared statistic", "p-value")
  colnames(t) <- ""
  
  res <- list(test.res = t,
              mu = object$mu,
              stdev = sqrt(object$var),
              bins = object$bins,
              df = object$parameter)
  
  class(res) <- "summary.ntestMeans"
  res
}

