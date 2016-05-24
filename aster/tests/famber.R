
 library(aster)

 ifam <- fam.bernoulli()

 p <- seq(0.1, 0.9, 0.1)
 theta <- log(p) - log(1 - p)

 zeroth <- double(length(p))
 first <- double(length(p))
 second <- double(length(p))

 for (i in seq(along = p)) {
    zeroth[i] <- famfun(ifam, 0, theta[i])
    first[i] <- famfun(ifam, 1, theta[i])
    second[i] <- famfun(ifam, 2, theta[i])
 }

 all.equal(zeroth, - log(1 - p))
 all.equal(first, p)
 all.equal(second, p * (1 - p))

