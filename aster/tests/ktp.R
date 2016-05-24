
 library(aster)

 do.chisq.test <- function(x, k, mu, max.bin) {
     stopifnot(all(x > k))
     stopifnot(k + 1 < max.bin)
     xx <- seq(k + 1, max.bin)
     yy <- dpois(xx, mu)
     yy[length(yy)] <- ppois(max.bin - 1, mu, lower.tail = FALSE)
     pp <- yy / sum(yy)
     ecc <- length(x) * pp
     if (any(ecc < 5.0))
         warning("violates rule of thumb about > 5 expected in each cell")
     cc <- tabulate(x, max.bin)
     cc <- cc[xx]
     cc[length(cc)] <- nsim - sum(cc[- length(cc)])
     chisqstat <- sum((cc - ecc)^2 / ecc)
     cat("chi squared statistic =", chisqstat, "\n")
     cat("degrees of freedom =", length(ecc) - 1, "\n")
     pval <- pchisq(chisqstat, length(ecc) - 1, lower.tail = FALSE)
     cat("p-value =", pval, "\n")
     if (exists("save.min.pval")) {
         save.min.pval <<- min(pval, save.min.pval)
         save.ntests <<- save.ntests + 1
     } else {
         save.min.pval <<- pval
         save.ntests <<- 1
     }
     foo <- rbind(cc, ecc)
     dimnames(foo) <- list(c("observed", "expected"), as.character(xx))
     print(foo)
 }

 set.seed(42)
 nsim <- 1e4

 mu <- 10
 k <- 2
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 22)

 mu <- 3.5
 k <- 2
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 11)

 mu <- 2.5
 k <- 2
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 10)

 mu <- 1.5
 k <- 2
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 8)

 mu <- 0.5
 k <- 2
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 6)

 nsim <- 1e5
 mu <- 0.1
 k <- 2
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 5)

 mu <- 0.01
 k <- 2
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 4)

 mu <- 1.5
 xpred <- 0:10
 save.seed <- .Random.seed
 x <- rktp(xpred, k, mu, xpred)
 .Random.seed <- save.seed
 my.x <- rep(0, length(xpred))
 for (i in seq(along = xpred))
     if (xpred[i] > 0)
         for (j in 1:xpred[i])
             my.x[i] <- my.x[i] + rktp(1, k, mu)
 all.equal(x, my.x)

 k <- 5
 mu <- pi
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 14)

 k <- 10
 mu <- exp(2)
 x <- rktp(nsim, k, mu)
 do.chisq.test(x, k, mu, 22)

 cat("number of tests:", save.ntests, "\n")
 cat("minimum p-value:", save.min.pval, "\n")
 cat("Bonferroni corrected minimum p-value:",
     save.ntests * save.min.pval, "\n")

 #####

 set.seed(42)
 nind <- 25
 nnode <- 1
 ncoef <- 1
 k <- 2

 pred <- 0
 fam <- 4

 theta <- 4 / 3
 mu <- exp(theta)
 x <- rpois(100, mu)
 x <- x[x > k]
 x <- x[1:nind]

 modmat <- matrix(1, nrow = nind, ncol = 1)
 out <- mlogl(theta, pred, fam, x, modmat, modmat, deriv = 2,
     type = "conditional")
 print(out)

 xxx <- seq(0, 100)
 ppp <- dpois(xxx, mu)
 ppp[xxx <= k] <- 0
 ppp <- ppp / sum(ppp)
 tau <- sum(xxx * ppp)

 sum(x - tau)
 max(abs( sum(x - tau) + out$gradient )) < 1e-12

 length(x) * sum((xxx - tau)^2 * ppp)

 epsilon <- 1e-6
 oute <- mlogl(theta + epsilon, pred, fam, x, modmat, modmat, deriv = 2,
     type = "conditional")
 (oute$value - out$value) / epsilon
 all.equal((oute$value - out$value) / epsilon, out$gradient,
     tol = 10 * epsilon)

 (oute$gradient - out$gradient) / epsilon
 all.equal((oute$gradient - out$gradient) / epsilon, as.numeric(out$hessian),
     tol = 10 * epsilon)

