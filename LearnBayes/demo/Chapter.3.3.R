##########################################################
# Section 3.3 Estimating a Heart Transplant Mortality Rate
##########################################################

 alpha=16;beta=15174
 yobs=1; ex=66
 y=0:10
 lam=alpha/beta
 py=dpois(y, lam*ex)*dgamma(lam, shape = alpha,
   rate = beta)/dgamma(lam, shape= alpha + y,
   rate = beta + ex)
 cbind(y, round(py, 3))

 lambdaA = rgamma(1000, shape = alpha + yobs, rate = beta + ex)

 ex = 1767; yobs=4
 y = 0:10
 py = dpois(y, lam * ex) * dgamma(lam, shape = alpha, 
     rate = beta)/dgamma(lam, shape = alpha + y,
     rate = beta + ex)
 cbind(y, round(py, 3))

 lambdaB = rgamma(1000, shape = alpha + yobs, rate = beta + ex)

 par(mfrow = c(2, 1))
 plot(density(lambdaA), main="HOSPITAL A", xlab="lambdaA", lwd=3)
 curve(dgamma(x, shape = alpha, rate = beta), add=TRUE)
 legend("topright",legend=c("prior","posterior"),lwd=c(1,3))
 plot(density(lambdaB), main="HOSPITAL B", xlab="lambdaB", lwd=3)
 curve(dgamma(x, shape = alpha, rate = beta), add=TRUE)
 legend("topright",legend=c("prior","posterior"),lwd=c(1,3))
