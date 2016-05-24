fitstatistics <- function (x,y,newx,newy,n,method,period){
if (method == "harmonic2") {
  d.f. <- 6
 sigmadf <- 2
}
else {
  d.f. <- 5
  sigmadf <- 1
}
MFPE <- crossprod(cbind(y,x) - cbind(newy,newx))/(n - d.f.)
Sigma <- MFPE*(n - d.f.)/n
AIC.y <- log(Sigma[1,1])*n + 4 * 2
AIC <- log(prod(diag(Sigma)))*n + (d.f.+sigmadf)*2
  #Based on multivariate AIC formula with covariance terms replaced.
return(c("MFPE" = MFPE, "AIC.y"=AIC.y,
         "AIC" = AIC,"n"=n, "period"=period,"d.f."=n-d.f.))
}
