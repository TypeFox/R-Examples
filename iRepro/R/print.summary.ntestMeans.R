print.summary.ntestMeans <-
function(x, ...){
  cat("Pearson's Chi-squared test for normality of means\n\n")
  cat("Baseline distribution: normal with mean",x$mu," and standard deviation",x$stdev,"\n\n")
  cat("Number of categories: ",x$bins,"\n")
  cat("Degrees of freedom: ",x$df,"\n")
  print(x$test.res)  
}

