
 library(aster)

 alpha <- 2.222
 ifam <- fam.negative.binomial(alpha)

 ##### usual size theta #####

 p <- seq(0.9, 0.1, -0.1)
 theta <- log(1 - p)
 qq <- exp(theta)
 pp <- (- expm1(theta))
 all.equal(p, pp)
 all.equal(pp, 1 - qq)

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, alpha * (- log(pp)))
 all.equal(first, alpha * qq / pp)
 all.equal(second, alpha * qq / pp^2)

 ##### very large negative theta #####

 rm(p)

 theta <- seq(-100, -10, 10)
 qq <- exp(theta)
 pp <- (- expm1(theta))

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, alpha * (- log(pp)))
 all.equal(first, alpha * qq / pp)
 all.equal(second, alpha * qq / pp^2)

 ##### very small negative theta #####

 theta <- (- 10^(- c(1:9, seq(10, 100, 10))))
 qq <- exp(theta)
 pp <- (- expm1(theta))

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, alpha * (- log(pp)))
 all.equal(first, alpha * qq / pp)
 all.equal(second, alpha * qq / pp^2)

 ##### random #####

 nind <- 50
 theta <- rep(- 1.75, nind)
 pred <- 0
 fam <- 1
 root <- seq(1, by = 0.5, length = nind)
 theta <- cbind(theta)
 root <- cbind(root)

 set.seed(42)
 rout <- raster(theta, pred, fam, root, famlist = list(ifam))

 set.seed(42)
 rout.too <- rnbinom(nind, size = alpha * as.numeric(root),
     prob = as.numeric(1 - exp(theta)))

 all.equal(as.numeric(rout), rout.too)



