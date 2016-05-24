
library(potts)

set.seed(42)

ncolor <- as.integer(4)
beta <- log(1 + sqrt(ncolor))
theta <- c(rep(0, ncolor), beta)

nrow <- 100
ncol <- 100
x <- matrix(1, nrow = nrow, ncol = ncol)
foo <- packPotts(x, ncolor)

out <- potts(foo, theta, nbatch = 10)
out$batch

out <- potts(foo, theta, nbatch = 10, boundary = "free")
out$batch

out <- potts(foo, theta, nbatch = 10, boundary = "condition")
out$batch

save.seed <- .Random.seed
out <- potts(foo, theta, nbatch = 20, blen = 10, nspac = 5)
niter <- out$nbatch * out$blen * out$nspac
.Random.seed <- save.seed
out.too <- potts(foo, theta, nbatch = niter)

my.batch <- out.too$batch[1:niter %% out$nspac == 0, ]
my.batch <- array(as.vector(my.batch),
    dim = c(out$blen, out$nbatch, ncolor + 1))
my.batch <- apply(my.batch, c(2, 3), mean)
all.equal(my.batch, out$batch)


