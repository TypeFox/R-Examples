gr0 <-
function(para, dat){ 
  beta0 <- para[1]; v <- para[2]; Delta.i <- dat$Delta.i; X.i <- dat$X.i
  n <- length(X.i)
  A.i <- 1/v + Delta.i
  B.i <- 1/v + exp(beta0)*X.i
  l.beta0 <- sum(Delta.i) - exp(beta0)* sum(A.i*X.i/B.i)
  l.v <-  (-sum(digamma(A.i)) + n*digamma(1/v) - n + n*log(v) + sum(log(B.i)) + sum(A.i/B.i))/v^2 
  -c(l.beta0, l.v)
}
