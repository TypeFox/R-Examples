# syn.nbc.r 
# Table 10.9: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
#     with assistance of: Andrew Robinson, University of Melbourne, Australia
# Synthetic NB-C estimated using optimiation
#
library(MASS)
nobs <- 50000
x2 <- runif(nobs)
x1 <- runif(nobs)
a <- 1.15            # value of alpha: 1.15               
xb <- 1.25*x1 + .1*x2 - 1.5   
mu <- 1/((exp(-xb)-1)*a)
p <- 1/(1+a*mu)
r <- 1/a
nbcy <- rnbinom(50000, size=r, prob = p)

nbc.reg.ml <- function(b.hat, X, y) {
 a.hat <- b.hat[1]
 xb.hat <- X %*% b.hat[-1]
 mu.hat <- 1 / ((exp(-xb.hat)-1)*a.hat)
 p.hat <- 1 / (1 + a.hat*mu.hat)
 r.hat <- 1 / a.hat
 sum(dnbinom(y,
             size = r.hat,
             prob = p.hat,
             log = TRUE))
}
nbcX <- cbind(1, x1, x2)
 p.0 <- c(a.hat = 1,
        b.0 = -2,
        b.1 = 1,
        b.2 = 1)

fit <- optim(p.0,        ## Maximize the JCLL
            nbc.reg.ml,
            X = nbcX,
            y = nbcy,
            control = list(fnscale = -1),
            hessian = TRUE
            )
fit$par                      ## ML estimates
sqrt(diag(solve(-fit$hessian)))   ## SEs

