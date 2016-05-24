rm(list=ls())
load("c://source/medpar.RData")
medpar$type <- factor(medpar$type)

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
# Create the design matrix
nbcX <- model.matrix(~ hmo + white + type, data = medpar)

# Starting points (discovered by trial and error!)
p.0 <- c(alpha = 0.5,
         cons = -1,
         hmo = 0,
         white = 0,
         type2 = 0,
         type3 = 0)

# Maximize the joint conditional LL
nbc.fit <- optim(p.0,             
                 nbc.reg.ml,
                 X = nbcX,
                 y = medpar$los,
                 control = list(
                   fnscale = -1,
                   maxit = 10000),
                 hessian = TRUE
                 )

# and obtain the parameter estimates and asymptotic SE's by 
(nbc.beta.hat <- nbc.fit$par)
(nbc.se.beta.hat <- sqrt(diag(solve(-nbc.fit$hessian))))
nbc.results <- data.frame(Estimate = nbc.beta.hat,
                          SE = nbc.se.beta.hat,
                          Z = nbc.beta.hat / nbc.se.beta.hat,
                          LCL = nbc.beta.hat - 1.96 * nbc.se.beta.hat,
                          UCL = nbc.beta.hat + 1.96 * nbc.se.beta.hat)
rownames(nbc.results) <- c("alpha", colnames(nbcX))
nbc.results <- nbc.results[c(2:nrow(nbc.results), 1),]
nbc.results

