
jllm <- function(y, mu, m, a) UseMethod("jllm")
linkFn <- function(mu, m, a) UseMethod("linkFn")
lPrime <- function(mu, m, a) UseMethod("lPrime")
unlink <- function(y, eta, m, a) UseMethod("unlink")
variance <- function(mu, m, a) UseMethod("variance")


devianceResids <- function(y, mu, m = 1, a = 1)
  sign(y - mu) * sqrt(2 * abs(jllm(y, mu, m, a) -
                              jllm(y, y,  m, a)))

devIRLS <- function(object, ...)
   sum(devianceResids(object, ...)^2)

initialize <- function(y, m = 1) {
   ret.y <- rep(mean(y), length(y))
   class(ret.y) <- class(y)
   ret.y
}

jllm.binomial <- function(y, mu, m = 1, a = 1) 
  dbinom(x = y, size = m, prob = mu / m, log = TRUE)
jllm.negBinomial <- function(y, mu, m = 1, a = 1) {
  dnbinom(y, mu = mu, size = 1 / a, log = TRUE)
}
jllm.poisson <- function(y, mu, m = 1, a = 1) {
  dpois(x = y, lambda = mu, log = TRUE)
}


variance.binomial <- function(mu, m = 1, a = 1) mu * (1 - mu/m)
variance.negBinomial <- function(mu, m = 1, a = 1) mu + a*mu^2
variance.poisson <- function(mu, m = 1, a = 1) mu


linkFn.logit <- function(mu, m = 1, a = 1) log(mu / (m - mu))
lPrime.logit <- function(mu, m = 1, a = 1) m / (mu * (m - mu))
unlink.logit <- function(y, eta, m = 1, a = 1) m / (1 + exp(-eta))

linkFn.probit <- function(mu, m = 1, a = 1) qnorm(mu / m)
lPrime.probit <- function(mu, m = 1, a = 1) 1 / (m * dnorm(qnorm(mu/m)))
unlink.probit <- function(y, eta, m = 1, a = 1) m * pnorm(eta)

linkFn.log <- function(mu, m = 1, a = 1) log(mu)
lPrime.log <- function(mu, m = 1, a = 1) 1/mu
unlink.log <- function(y, eta, m = 1, a = 1) exp(eta)  

linkFn.identity <- function(mu, m = 1, a = 1) mu
lPrime.identity <- function(mu, m = 1, a = 1) 1
unlink.identity <- function(y, eta, m = 1, a = 1) eta

linkFn.negbin <- function(mu, m=1, a=1) -log(1+1/(a*mu))
lPrime.negbin <- function(mu, m=1, a=1) 1/(mu+a*mu^2)
unlink.negbin <- function(y, eta, m=1, a=1) 1/(a*(exp(-eta) - 1))

linkFn.cloglog <- function(mu, m=1, a=1) log(-log(m - mu))
lPrime.cloglog <- function(mu, m=1, a=1) 1/((mu - m)*log(1 - mu/m))
unlink.cloglog <- function(y, eta, m=1, a=1) m - exp(-exp(-eta))

linkFn.cloglog <- function(mu, m=1, a=1) log(-log(1-mu/m))
lPrime.cloglog <- function(mu, m=1, a=1) 1/((mu - m)*log(1-mu/m))
unlink.cloglog <- function(y, eta, m=1, a=1) m*(1-exp(-exp(eta)))

linkFn.inverse <- function(mu, m=1, a=1) 1/mu
lPrime.inverse <- function(mu, m=1, a=1) -1/mu^2
unlink.inverse <- function(y, eta, m=1, a=1) 1/eta

linkFn.inverse2 <- function(mu, m=1, a=1) 1/mu^2
lPrime.inverse2 <- function(mu, m=1, a=1) -1/mu^3
unlink.inverse2 <- function(y, eta, m=1, a=1) 1/sqrt(eta)




#### Chapter 5

jll <- function(y, y.hat, ...) UseMethod("jll")

jll.bernoulli <- function(y, y.hat, ...) {
  dbinom(x = y, size = 1, prob = y.hat, log = TRUE)
}

Sjll <- function(b.hat, X, y, offset = 0, ...) {
  y.hat <- predict(y, b.hat, X, offset, ...)
  sum(jll(y, y.hat, ...))
}

## Changed arguments for logit1 

# unlink <- function(y, eta) UseMethod("unlink")
unlink.logit1 <- function(y, eta, m=1, a=1) 1 / (1 + exp(-eta))
unlink.cloglog1 <- function(y, eta, m=1, a=1) 1-exp(-exp(eta))
unlink.probit1 <- function(y, eta, m=1, a=1) pnorm(eta)

unlink.logit1 <- function(y, eta, ...) 1 / (1 + exp(-eta))
unlink.cloglog1 <- function(y, eta, ...) 1-exp(-exp(eta))
unlink.probit1 <- function(y, eta, ...) pnorm(eta)

maximize <- function(start, f, X, y, offset = 0, ...) {
  optim(par = start,           
        fn = f,
        X = X,
        y = y,
        offset = offset,
        method = "BFGS",
        control = list(
          reltol = 1e-16,
          fnscale = -1,
          maxit = 10000),
        hessian = TRUE,
        ...
        )
}

kickStart <- function(y, X, offset)
                             UseMethod("kickStart")

kickStart.default <- function(y, X, offset = 0) {
  coef(lm(I(y - offset) ~ X - 1))
}

jll.poisson <- function(y, y.hat, ...) {
  dpois(y, lambda = y.hat, log = TRUE)
}

unlink.log <- function(y, eta, m=1, a=1) exp(eta)  

kickStart.log <- function(y, X, offset = 0) {
  coef(lm(I((log(y + 0.1) - offset) ~ X - 1)))
}

jll.ztp <- function(y, y.hat, ...) 
  dpois(y, lambda = y.hat, log = TRUE) - log(1 - exp(-y.hat))

jll2 <- function(y, y.hat, scale, ...) UseMethod("jll2")

jll2.negBinomial <- function(y, y.hat, scale, ...) {
  dnbinom(y,
          mu = y.hat,
          size = 1 / scale,
          log = TRUE)
}

jll2.nb2 <- function(y, y.hat, scale) {
  dnbinom(y,
          size = scale,
          mu = y.hat, log = TRUE)
}


unlink_s <- function(y, eta) UseMethod("unlink_s")

getDispersion <- function(y, scale) UseMethod("getDispersion")
getDispersion.negBinomial <- function(y, scale) 1

getDispersion.nb2 <- function(y, scale) 1


#kickStart.log <- function(y, X, family, offset = NULL) {
#  coef(lm(log(y + 0.01) ~ X - 1, offset = offset))
#}

unlink_s.log_s <- function(y, eta) exp(eta)
unlink_s.identity_s <- function(y, eta) eta
unlink_s.inverse_s <- function(y, eta) 1 / eta

jll2.normal <- function(y, y.hat, scale, ...) {
  dnorm(y,
        mean = y.hat,
        sd = scale, log = TRUE)
}
getDispersion.normal <- function(y, scale) scale^2

# unlink.identity <- function(y, eta) eta


jll2.gamma <- function(y, y.hat, scale, ...) {
  dgamma(y,
         shape = 1 / scale,
         scale = y.hat * scale, log = TRUE)
}

getDispersion.gamma <- function(y, scale) scale


pearsonResiduals2 <- function(y, b.hat, X, p, offset = 0) {
  y.hat <- predict(y, b.hat[1:p], X[,1:p], offset)
  scale <- predict_s(y, b.hat[-(1:p)], X[,-(1:p)])
  (y - y.hat) / sqrt(y.hat+ 1/scale*y.hat*y.hat) 
}


kickStart.inverse <- function(y, X, family, offset = NULL) {
  coef(lm(I(1/(y)) ~ X - 1), offset = offset)
}

