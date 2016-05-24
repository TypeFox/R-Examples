
 library(aster)

 do.chisq.test <- function(x, alpha, k, mu, max.bin) {
     stopifnot(all(x > k))
     stopifnot(k + 1 < max.bin)
     xx <- seq(k + 1, max.bin)
     yy <- dnbinom(xx, size = alpha, mu = mu)
     yy[length(yy)] <- pnbinom(max.bin - 1, size = alpha, mu = mu,
         lower.tail = FALSE)
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

 alpha <- 2.222
 mu <- 10
 k <- 2
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 40)

 alpha <- 2.222
 mu <- 3.5
 k <- 2
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 20)

 alpha <- 2.222
 mu <- 2.5
 k <- 2
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 16)

 alpha <- 2.222
 mu <- 1.5
 k <- 2
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 12)

 alpha <- 2.222
 mu <- 0.5
 k <- 2
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 8)

 alpha <- 2.222
 mu <- 0.1
 k <- 2
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 5)

 nsim <- 2e5
 alpha <- 2.222
 mu <- 0.01
 k <- 2
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 5)

 alpha <- 2.222
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

 nsim <- 1e4
 alpha <- 5.55
 k <- 5
 mu <- pi
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 16)

 alpha <- 5.55
 k <- 10
 mu <- exp(2)
 x <- rktnb(nsim, alpha, k, mu)
 do.chisq.test(x, alpha, k, mu, 29)

 cat("number of tests:", save.ntests, "\n")
 cat("minimum p-value:", save.min.pval, "\n")
 cat("Bonferroni corrected minimum p-value:",
     save.ntests * save.min.pval, "\n")

 #####

 set.seed(42)
 nind <- 25
 nnode <- 1
 ncoef <- 1
 alpha <- 3.333
 k <- 2

 pred <- 0
 fam <- 1
 ifam <- fam.truncated.negative.binomial(size = alpha, trunc = k)
 aster:::setfam(list(ifam))
 theta.origin <- aster:::getfam()[[1]]$origin

 theta <- (- 4 / 3)
 p <- 1 - exp(theta)
 x <- rnbinom(1000, size = alpha, prob = p)
 x <- x[x > k]
 x <- x[1:nind]
 modmat <- matrix(1, nrow = nind, ncol = 1)

 out <- mlogl(theta - theta.origin, pred, fam, x, modmat, modmat,
     deriv = 2, type = "conditional", famlist = list(ifam))
 print(out)

 xxx <- seq(0, 100)
 ppp <- dnbinom(xxx, size = alpha, prob = p)
 ppp[xxx <= k] <- 0
 ppp <- ppp / sum(ppp)
 tau <- sum(xxx * ppp)

 my.grad.logl <- sum(x - tau)
 all.equal(- out$gradient, my.grad.logl)

 my.fish.info <- length(x) * sum((xxx - tau)^2 * ppp)
 all.equal(as.numeric(out$hessian), my.fish.info)

