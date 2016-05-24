
qloggamma <- function(p,lambda,zero=0.001){
if (lambda < -zero) p <- 1-p
if (abs(lambda) > zero) {
  k   <- 1/lambda^2
  res <- (log(qgamma(p,shape=k,rate=1))+log(lambda^2))/lambda} 
else {res <- qnorm(p,mean=0,sd=1)} 
res}
 
