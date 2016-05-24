### This example is analog to "boot_par.r", and one can run it by the command
### SHELL> mpiexec -np 2 Rscript --vanilla 03_boot_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.rank <- comm.rank()
set.seed(123 + .comm.rank + 1)

library(boot)
cd4.rg <- function(data, mle) MASS::mvrnorm(nrow(data), mle$m, mle$v)
cd4.mle <- list(m = colMeans(cd4), v = var(cd4))
res <- boot(cd4, corr, R = 500, sim = "parametric",
            ran.gen = cd4.rg, mle = cd4.mle)

cd4.boot <- do.call(c, allgather(res))
ci <- boot.ci(cd4.boot,  type = c("norm", "basic", "perc"),
              conf = 0.9, h = atanh, hinv = tanh)
comm.print(ci)

finalize()
