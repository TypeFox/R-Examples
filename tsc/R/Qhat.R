Qhat <-
function(teststats, alpha, nsims, u0,sigma0){
  Jmin = floor(max(c(1,((1-alpha)*nsims-sqrt(nsims)*log(nsims)))))
  Jmax = ceiling(min(c(nsims,((1-alpha)*nsims+sqrt(nsims)*log(nsims)))))
  Jint = as.numeric(Jmin:Jmax)
  num = exp((-nsims/(2*alpha*(1-alpha)))*((1-alpha - Jint/nsims)^2))*(sqrt(sigma0/(2*pi))*(exp(-((teststats[Jint-1]-u0)^2)/(2*sigma0))-exp(-((teststats[Jint]-u0)^2)/(2*sigma0)))+u0*(pnorm(teststats[Jint], u0, sqrt(sigma0))-pnorm(teststats[Jint-1], u0, sqrt(sigma0))))
  den = exp((-nsims/(2*alpha*(1-alpha)))*((1-alpha - Jint/nsims)^2))*((pnorm(teststats[Jint], u0, sqrt(sigma0))-pnorm(teststats[Jint-1], u0, sqrt(sigma0))))
  value.cal = sum(num)/sum(den)
  return(value.cal)
}
