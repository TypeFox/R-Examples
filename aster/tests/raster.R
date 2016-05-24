
 library(aster)

 set.seed(42)
 nind <- 25
 nnode <- 5

 theta <- rnorm(nind * nnode) / 10
 theta <- matrix(theta, nind, nnode)

 root <- sample(1:3, nind * nnode, replace = TRUE)
 root <- matrix(root, nind, nnode)

 fam <- c(1, 1, 2, 3, 3)
 pred <- c(0, 1, 1, 2, 3)
 famnam <- sapply(fam.default(), as.character)

 save.seed <- .Random.seed
 x <- raster(theta, pred, fam, root)
 
 .Random.seed <- save.seed
 my.x <- NaN * x

 for (i in 1:nnode) {
    ipred <- pred[i]
    if (ipred == 0) {
        xpred <- root[ , i]
    } else {
        xpred <- my.x[ , ipred]
    }
    if (famnam[fam[i]] == "bernoulli") {
        p <- 1 / (1 + exp(- theta[ , i]))
        xi <- rep(0, nind)
        ii <- xpred > 0
        xi[ii] <- rbinom(sum(ii), xpred[ii], p[ii])
        my.x[ , i] <- xi
    }
    if (famnam[fam[i]] == "poisson") {
        mu <- exp(theta[ , i])
        xi <- rep(0, nind)
        ii <- xpred > 0
        xi[ii] <- rpois(sum(ii), xpred[ii] * mu[ii])
        my.x[ , i] <- xi
    }
    if (famnam[fam[i]] == "truncated.poisson(truncation = 0)") {
        mu <- exp(theta[ , i])
        my.x[ , i] <- rnzp(nind, mu, xpred)
    }
 }

 all.equal(x, my.x)

