
 library(aster)

 do.chisq.test <- function(x, mu, max.bin) {
     stopifnot(all(x >= 0))
     xx <- seq(1, max.bin)
     yy <- dpois(xx, mu)
     yy[length(yy)] <- ppois(max.bin - 1, mu, lower.tail = FALSE)
     pp <- yy / sum(yy)
     ecc <- length(x) * pp
     if (any(ecc < 5.0))
         warning("violates rule of thumb about > 5 expected in each cell")
     cc <- tabulate(x, max.bin)
     chisqstat <- sum((cc - ecc)^2 / ecc)
     cat("chi squared statistic =", chisqstat, "\n")
     cat("degrees of freedom =", length(ecc) - 1, "\n")
     cat("p-value =", pchisq(chisqstat, length(ecc) - 1, lower.tail = FALSE),
         "\n")
     foo <- rbind(cc, ecc)
     dimnames(foo) <- list(c("observed", "expected"), as.character(xx))
     print(foo)
 }

 set.seed(42)
 nsim <- 1e4

 mu <- 2.0
 x <- rnzp(nsim, mu)
 do.chisq.test(x, mu, 8)

 mu <- 1.0
 x <- rnzp(nsim, mu)
 do.chisq.test(x, mu, 5)

 mu <- 0.5
 x <- rnzp(nsim, mu)
 do.chisq.test(x, mu, 4)

 nsim <- 1e6
 mu <- 0.05
 x <- rnzp(nsim, mu)
 do.chisq.test(x, mu, 4)

 # nsim <- 1e7
 # mu <- 0.005
 # x <- rnzp(nsim, mu)
 # do.chisq.test(x, mu, 3)

 mu <- 0.5
 xpred <- 0:10
 save.seed <- .Random.seed
 x <- rnzp(xpred, mu, xpred)
 .Random.seed <- save.seed
 my.x <- rep(0, length(xpred))
 for (i in seq(along = xpred))
     if (xpred[i] > 0)
         for (j in 1:xpred[i])
             my.x[i] <- my.x[i] + rnzp(1, mu)
 all.equal(x, my.x)


