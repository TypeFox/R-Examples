# cramer : cramer.test

CvM.test <- function(x,y,...) {
  dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
  test <- cramer::cramer.test(x,y,...)
  names(test$statistic) <- "T"
  result <- list(statistic=test$statistic,p.value=test$p.value,alternative="two.sided",
    method="Two-sample Cram\u00E9r-von Mises test",data.name=dname)
  class(result) <- "htest"
  return(result)
}
