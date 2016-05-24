### This is an example modified from the help page of parLapply()
### One can source this file into R interactive model or run it by the command
### SHELL> Rscript --vanilla 03_boot_spmd.r

library(parallel)
cl <- makeCluster(mc <- getOption("cl.cores", 2))

run1 <- function(i) {
  set.seed(123 + i)

  library(boot)
  cd4.rg <- function(data, mle) MASS::mvrnorm(nrow(data), mle$m, mle$v)
  cd4.mle <- list(m = colMeans(cd4), v = var(cd4))
  boot(cd4, corr, R = 500, sim = "parametric",
       ran.gen = cd4.rg, mle = cd4.mle)
}

library(boot)
cd4.boot <- do.call(c, parLapply(cl, seq_len(mc), run1))
boot.ci(cd4.boot,  type = c("norm", "basic", "perc"),
        conf = 0.9, h = atanh, hinv = tanh)

stopCluster(cl)
