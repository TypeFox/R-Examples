
 library(aster)

 sigma <- 2.222
 ifam <- fam.normal.location(sigma)

 ##### usual size theta #####

 theta <- seq(-3.0, 3.0, 0.1)

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, sigma^2 * theta^2 / 2)
 all.equal(first, sigma^2 * theta)
 all.equal(second, sigma^2 * theta^0)

 ##### random #####

 theta <- seq(-3.5, 3.5, 0.1)
 nind <- length(theta)
 pred <- 0
 fam <- 1
 root <- seq(1, 5, length = nind)
 theta <- cbind(theta)
 root <- cbind(root)

 set.seed(42)
 rout <- raster(theta, pred, fam, root, famlist = list(ifam))

 set.seed(42)
 moo <- sigma^2 * theta * as.numeric(root)
 cow <- sigma * sqrt(as.numeric(root))
 rout.too <- rnorm(nind, mean = moo, sd = cow)

 all.equal(as.numeric(rout), rout.too)

