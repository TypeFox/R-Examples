### This example is analog to "glm_par.r", and one can run it by the command
### SHELL> mpiexec -np 2 Rscript --vanilla 04_glm_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()

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

time.proc <- system.time({
  id <- get.jid(length(obs))
  ret <- allgather(sapply(id, fitmodel, obs), unlist = TRUE)
})
comm.print(time.proc)

finalize()
