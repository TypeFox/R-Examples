plot.interv_detect <- function(x, ...){
  main <- "Test statistic over time"
  timser <- x$fit_H0$ts
  if(is.ts(timser)){
    plot(x=time(timser)[as.numeric(names(x$test_statistic_tau))], y=x$test_statistic_tau, main=main, xlab="Time", ylab="Test statistic", type ="o", ...)
  }else{
    plot(x=as.numeric(names(x$test_statistic_tau)), y=x$test_statistic_tau, main=main, xlab="Time", ylab="Test statistic", type="o", ...)
  }
  invisible()
}
