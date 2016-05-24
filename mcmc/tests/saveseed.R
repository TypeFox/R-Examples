
 library(mcmc)

 set.seed(42)

 h <- function(x) if (all(x >= 0) && sum(x) <= 1) return(1) else return(-Inf)
 out <- metrop(h, initial = rep(0, 5), nbatch = 100, blen = 17, nspac = 3,
     scale = 0.1)

 save.seed <- .Random.seed

 out1 <- metrop(out)
 out2 <- metrop(out1)
 out3 <- metrop(out, nbatch = 2 * out$nbatch)

 fred <- rbind(out1$batch, out2$batch)
 identical(fred, out3$batch)

