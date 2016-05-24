# p.reg.ml.r
# Poisson optimization  
# Ch 6.2.1 : Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
set.seed(3357)
b <- c(5, 1, 0.5)                    ## Population parameters
n <- 10000
X <- cbind(1, rnorm(n), rnorm(n))    ## Design matrix
y <- rpois(n = n, lambda = X %*% b)

p.reg.ml <- function(b.hat, X, y) {  ## Joint Conditional LL
  sum(dpois(y, lambda = X %*% b.hat, log = TRUE))  
}
p.0 <- lm.fit(X, y)$coef             ## Obtain initial estimates
fit <- optim(p.0,                    ## Maximize JCLL
             p.reg.ml,
             X = X,
             y = y,
             control = list(fnscale = -1),
             hessian = TRUE
             )
stderr <- sqrt(diag(solve(-fit$hessian))) ## Asymptotic SEs
poiresults <- data.frame(fit$par, stderr)
poiresults
