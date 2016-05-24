
 library(aster)

 ifam <- fam.truncated.poisson(truncation = 0)

 mu <- seq(0.1, 3.0, 0.2)
 theta <- log(mu)

 zeroth <- double(length(theta))
 first <- double(length(theta))
 second <- double(length(theta))

 for (i in seq(along = theta)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, log(exp(mu) - 1))
 tau <- mu / (1 - exp(- mu))
 all.equal(first, tau)
 all.equal(second, tau * (1 - tau * exp(- mu)))

