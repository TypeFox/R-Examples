
 library(aster)

 alpha <- 2.222
 k <- 2
 ifam <- fam.truncated.negative.binomial(alpha, k)
 aster:::setfam(list(ifam))
 aster:::getfam()
 aster:::clearfam()

 ##### usual size theta #####

 p <- seq(0.9, 0.1, -0.1)
 theta <- log(1 - p)
 qq <- exp(theta)
 pp <- (- expm1(theta))
 all.equal(p, pp)
 all.equal(pp, 1 - qq)
 beeta <- pnbinom(k + 1, size = alpha, prob = pp, lower.tail = FALSE) /
     dnbinom(k + 1, size = alpha, prob = pp)

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, alpha * (- log(pp)) +
     pnbinom(k, size = alpha, prob = pp, lower.tail = FALSE, log.p = TRUE))
 all.equal(first, alpha * qq / pp + (k + 1) / (1 + beeta) / pp)
 all.equal(second, alpha * qq / pp^2 - (k + 1) / (1 + beeta) / pp^2 *
     (- qq + (k + 1 + alpha) * qq / (1 + beeta) +
     (alpha - pp * (k + 1 + alpha)) * beeta / (1 + beeta)))

 ##### usual size theta (continued) #####
 ##### check by numerical derivative #####

 epsilon <- 1e-8

 zeroth.minus <- zeroth
 zeroth.plus <- zeroth
 for (i in seq(along = theta)) {
     zeroth.minus[i] <- famfun(ifam, 0, theta[i] - epsilon)
     zeroth.plus[i] <- famfun(ifam, 0, theta[i] + epsilon)
 }
 all.equal(first, (zeroth.plus - zeroth.minus) / (2 * epsilon),
     tolerance = sqrt(epsilon))

 first.minus <- first
 first.plus <- first
 for (i in seq(along = theta)) {
     first.minus[i] <- famfun(ifam, 1, theta[i] - epsilon)
     first.plus[i] <- famfun(ifam, 1, theta[i] + epsilon)
 }
 all.equal(second, (first.plus - first.minus) / (2 * epsilon),
     tolerance = sqrt(epsilon))

 ##### very large negative theta #####

 rm(p)

 theta <- seq(-100, -10, 10)
 qq <- exp(theta)
 pp <- (- expm1(theta))
 beeta.up <- pnbinom(k + 1, size = alpha, prob = pp, lower.tail = FALSE)
 beeta.dn <- dnbinom(k + 1, size = alpha, prob = pp)
 beeta <- beeta.up / beeta.dn
 beeta[beeta.up == 0] <- 0

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, alpha * (- log(pp)) +
     pnbinom(k, size = alpha, prob = pp, lower.tail = FALSE, log.p = TRUE))
 all.equal(first, alpha * qq / pp + (k + 1) / (1 + beeta) / pp)
 all.equal(second, alpha * qq / pp^2 - (k + 1) / (1 + beeta) / pp^2 *
     (- qq + (k + 1 + alpha) * qq / (1 + beeta) +
     (alpha - pp * (k + 1 + alpha)) * beeta / (1 + beeta)))

 ##### very small negative theta #####

 theta <- (- 10^(- c(1:9, seq(10, 100, 10))))
 qq <- exp(theta)
 pp <- (- expm1(theta))
 beeta <- pnbinom(k + 1, size = alpha, prob = pp, lower.tail = FALSE) /
     dnbinom(k + 1, size = alpha, prob = pp)

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, alpha * (- log(pp)) +
     pnbinom(k, size = alpha, prob = pp, lower.tail = FALSE, log.p = TRUE))
 all.equal(first, alpha * qq / pp + (k + 1) / (1 + beeta) / pp)
 all.equal(second, alpha * qq / pp^2 - (k + 1) / (1 + beeta) / pp^2 *
     (- qq + (k + 1 + alpha) * qq / (1 + beeta) +
     (alpha - pp * (k + 1 + alpha)) * beeta / (1 + beeta)))

 ##### random #####

 nind <- 50
 theta <- rep(- 1.75, nind)
 pred <- 0 
 fam <- 1
 root <- rep(1:3, length.out = nind)
 theta <- cbind(theta)
 root <- cbind(root)
 qq <- exp(theta)
 pp <- (- expm1(theta))
 mu <- alpha * qq / pp
    
 set.seed(42)
 rout <- raster(theta, pred, fam, root, famlist = list(ifam))
 
 set.seed(42)
 rout.too <- rktnb(nind, alpha, k, mu, root)
 
 all.equal(as.numeric(rout), rout.too)

