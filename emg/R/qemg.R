qemg <- function(p, mu=0, sigma=1, lambda=1, lower.tail=TRUE, log.p=FALSE)
{
  if(length(p[p>1]) > 1) {stop("p must be equal to or less than 1")}
  if(length(p[p<0]) > 1) {stop("p must be equal to or greater than 0")}
  
  if(!lower.tail) {p <- 1-p}
  
  result <- rep(NaN, length(p))
  for(i in 1:length(p))
  {
    result[i] <- optim(c(mu),
          fn=function(y) {abs(p[i] - pemg(y, mu, sigma, lambda))},
          lower=-Inf, upper=Inf, method="L-BFGS-B")$par
  }
  
  if(log.p) {log(result)} else {result}
}
