# nb.reg.ml.r  
# Synthetic MLE Negative Binomial NB2 data and model 
# Table 9.11: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
#     with assistance of: Andrew Robinson, University of Melbourne, Australia
#
set.seed(85132)
b <- c(5, 2, 3, 0.5)             ## Population parameters
n <- 10000
X <- cbind(rlnorm(n), rlnorm(n)) ## Design matrix
y <- rnbinom(n = n,              ## Choice of parameterization
             mu = b[1] + b[2] * X[,1],
             size = b[3] + b[4] * X[,2])
nb.reg.ml <- function(b.hat, X, y) {  ## JCLL
  sum(dnbinom(y,
              mu = b.hat[1] + b.hat[2] * X[,1],
              size = b.hat[3] + b.hat[4] * X[,2],
              log = TRUE))  
}
p.0 <- c(1,1,1,1)                ## initial estimates
fit <- optim(p.0,                ## Maximize the JCLL
             nb.reg.ml,
             X = X,
             y = y,
             control = list(fnscale = -1),
             hessian = TRUE
             )
stderr <- sqrt(diag(solve(-fit$hessian))) ## Asymptotic SEs
nbresults <- data.frame(fit$par, stderr)
nbresults


















