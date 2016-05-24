plot.tsglm <- function(x, ask=TRUE, ...){
  devAskNewPage(ask=ask)
  residu <- residuals(x, type="pearson")
  acf(residu, main="ACF of Pearson residuals")
  #hist(residu, main="Histogram of Pearson residuals", xlab="Residuals")
  plot(residu, type="o", xlab="Time", ylab="Residuals", main="Pearson residuals over time")
  par(ann=FALSE)
  cpgram(residu, main="")
  title(main="Cumulative periodogram of Pearson residuals", xlab="Frequency")
  par(ann=TRUE)
  pit(x)
  marcal(x)
  invisible()
}
