BootC <-
function (Data) {
BSs <- function(bootdata, indices) {
  d <- bootdata[indices,] # allows boot to select sample 
  fit <- lm(y ~ x1 + x2 + x1*x2, data=d)
  coefficient<-summary(fit)$coefficients
  covariance<-vcov(fit)
  b2<-coefficient[3]
  b3<-coefficient[4]
  v22<-covariance[3,3]
  v23<-covariance[3,4]
  v33<-covariance[4,4]
  C <-(-1)*b2/b3
  results<-c(C)
  return(results)
} 
temp <- boot(data=Data, statistic=BSs, R=1000)
return(boot.ci(temp, type=c("norm","basic","perc","bca"), index=1))
}
