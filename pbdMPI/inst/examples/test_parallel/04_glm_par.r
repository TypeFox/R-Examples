### This is an example to apply glm() in a large number of tasks.
### One can source this file into R interactive model or run it by the command
### SHELL> Rscript --vanilla 04_glm_spmd.r

library(parallel)

rdata <- function(n){
  x <- 2 * runif(n) - 1
  m <- 0 + 1 * x + 2 * x^2 + 3 * x^3
  y <- exp(m)
  z <- round(exp(m + rnorm(n)), 0)
  data.frame(list(y = y, z = z, x1 = x, x2 = x^2, x3 = x^3))
}
fitmodel <- function(i, dat){
  mod <- glm(z ~ x1 + x2 + x3, data = dat[[i]], family = poisson)
  mod$coefficients
}
set.seed(123)
obs <- lapply(rep(2000, 500), rdata)

system.time(ret <- lapply(1:length(obs), fitmodel, obs))
system.time(ret.mc <- mclapply(1:length(obs), fitmodel, obs))
