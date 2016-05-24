summary.ntestRes <-
function(object, ...){
  t <- as.data.frame(cbind(c(object$statistic.res1,object$p.value.res1),c(object$statistic.res2,object$p.value.res2)))
  rownames(t) <- c("Chi-squared statistic", "p-value")
  colnames(t) <- c("time = 1", "time = 2")
  
  res <- list(test.res = t,
              mu = object$mu,
              stdev = sqrt(object$var),
              bins = object$bins,
              df = object$parameter)
  
  class(res) <- "summary.ntestRes"
  res
}

