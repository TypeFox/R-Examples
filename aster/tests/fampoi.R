
 library(aster)

 ifam <- fam.poisson()

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

 all.equal(zeroth, mu)
 all.equal(first, mu)
 all.equal(second, mu)

