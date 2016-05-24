
 library(aster2)

 set.seed(42)

 epsilon <- 1e-6

 ##### Bernoulli #####

 theta <- seq(-3, 3, 0.1)

 cumfun <- function(theta)
     ifelse(theta > 0, theta + log1p(exp(- theta)), log1p(exp(theta)))

 moofun <- function(theta)
     ifelse(theta > 0, 1 / (1 + exp(- theta)), exp(theta) / (1 + exp(theta)))

 voofun <- function(theta) {
     pp <- moofun(theta)
     qq <- moofun(- theta)
     pp * qq
 }

 thoofun <- function(theta) {
     pp <- moofun(theta)
     qq <- moofun(- theta)
     - pp * qq * tanh(theta / 2)
 }

 zeroth <- mapply(function(theta) cumulant(theta, fam.bernoulli())$zeroth,
     theta)

 all.equal(zeroth, cumfun(theta))

 first <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     deriv = 1)$first, theta)

 all.equal(first, moofun(theta))

 second <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     deriv = 2)$second, theta)

 all.equal(second, voofun(theta))

 third <- mapply(function(theta) cumulant(theta, fam.bernoulli(),
     deriv = 3)$third, theta)

 all.equal(third, thoofun(theta))

 grad <- moofun(theta)
 fdgrad <- (cumfun(theta + epsilon) - cumfun(theta - epsilon)) / (2 * epsilon)
 all.equal(grad, fdgrad)

 hess <- voofun(theta)
 fdhess <- (moofun(theta + epsilon) - moofun(theta - epsilon)) / (2 * epsilon)
 all.equal(hess, fdhess)

 thss <- thoofun(theta)
 fdth <- (voofun(theta + epsilon) - voofun(theta - epsilon)) / (2 * epsilon)
 all.equal(thss, fdth)

 foo <- mapply(function(xi) link(xi, fam.bernoulli(), deriv = 0)$zeroth, first)

 all.equal(theta, foo)

 foo <- mapply(function(xi) link(xi, fam.bernoulli(), deriv = 1)$first, first)

 all.equal(second, 1 / foo)

 ##### Poisson #####

 zeroth <- mapply(function(theta) cumulant(theta, fam.poisson())$zeroth, theta)

 all.equal(zeroth, exp(theta))

 first <- mapply(function(theta) cumulant(theta, fam.poisson(),
     deriv = 1)$first, theta)

 all.equal(first, exp(theta))

 second <- mapply(function(theta) cumulant(theta, fam.poisson(),
     deriv = 2)$second, theta)

 all.equal(second, exp(theta))

 third <- mapply(function(theta) cumulant(theta, fam.poisson(),
     deriv = 3)$third, theta)

 all.equal(third, exp(theta))

 foo <- mapply(function(xi) link(xi, fam.poisson(), deriv = 0)$zeroth, first)

 all.equal(theta, foo)

 foo <- mapply(function(xi) link(xi, fam.poisson(), deriv = 1)$first, first)

 all.equal(second, 1 / foo)

 ##### zero-truncated Poisson #####

 cumfun <- function(theta) {
     m <- exp(theta)
     result <- m + log1p(- exp(- m))
     bar <- m / 2 * (1 + m / 3 * (1 + m / 4 * (1 + m / 5 * (1 + m / 6 *
         (1 + m / 7 * (1 + m / 8))))))
     inies <- theta < (- 4)
     result[inies] <- (theta + log1p(bar))[inies]
     return(result)
 }

 moofun <- function(theta) {
     m <- exp(theta)
     result <- m / (1 - exp(- m))
     bar <- m / 2 * (1 + m / 3 * (1 + m / 4 * (1 + m / 5 * (1 + m / 6 *
         (1 + m / 7 * (1 + m / 8))))))
     inies <- theta < (- 4)
     result[inies] <- (m + 1 / (1 + bar))[inies]
     return(result)
 }

 voofun <- function(theta) {
     m <- exp(theta)
     mu <- moofun(theta)
     return(mu * (1 + m - mu))
 }

 thoofun <- function(theta) {
     m <- exp(theta)
     mu <- moofun(theta)
     return(mu * ((1 + m - mu) * (1 + m - 2 * mu) + m))
 }

 theta <- seq(-12, 2, 0.1)

 zeroth <- mapply(function(theta)
     cumulant(theta, fam.zero.truncated.poisson())$zeroth, theta)

 all.equal(zeroth, cumfun(theta))

 first <- mapply(function(theta) cumulant(theta, fam.zero.truncated.poisson(),
     deriv = 1)$first, theta)

 all.equal(first, moofun(theta))

 second <- mapply(function(theta) cumulant(theta, fam.zero.truncated.poisson(),
     deriv = 2)$second, theta)

 all.equal(second, voofun(theta))

 third <- mapply(function(theta) cumulant(theta, fam.zero.truncated.poisson(),
     deriv = 3)$third, theta)

 all.equal(third, thoofun(theta))

 grad <- moofun(theta)
 fdgrad <- (cumfun(theta + epsilon) - cumfun(theta - epsilon)) / (2 * epsilon)
 all.equal(grad, fdgrad)

 hess <- voofun(theta)
 fdhess <- (moofun(theta + epsilon) - moofun(theta - epsilon)) / (2 * epsilon)
 all.equal(hess, fdhess)

 thss <- thoofun(theta)
 fdth <- (voofun(theta + epsilon) - voofun(theta - epsilon)) / (2 * epsilon)
 all.equal(thss, fdth)

 foo <- mapply(function(xi) link(xi, fam.zero.truncated.poisson(),
     deriv = 0)$zeroth, first)

 all.equal(theta, foo)

 foo <- mapply(function(xi) link(xi, fam.zero.truncated.poisson(),
     deriv = 1)$first, first)

 all.equal(second, 1 / foo)

 ##### Multinomial #####

 d <- 4

 cumfun <- function(theta) {
     stopifnot(length(theta) == d)
     theta.max = max(theta)
     j <- seq(along = theta)[theta == theta.max]
     j <- j[1]
     return(theta.max + log1p(sum(exp(theta - theta.max)[-j])))
 }

 moofun <- function(theta) {
     stopifnot(length(theta) == d)
     theta.max = max(theta)
     return(exp(theta - theta.max) / sum(exp(theta - theta.max)))
 }

 voofun <- function(theta) {
     mu <- moofun(theta)
     diag(mu) - outer(mu, mu)
 }

 thoofun <- function(theta) {
     mu <- moofun(theta)
     result <- 2 * outer(mu, outer(mu, mu))
     for (i in 1:d) {
         for (j in 1:d) {
             if (i != j) {
                 qux <- mu[i] * mu[j]
                 result[i, i, j] <- result[i, i, j] - qux
                 result[i, j, i] <- result[i, j, i] - qux
                 result[j, i, i] <- result[j, i, i] - qux
             }
         }
         qux <- mu[i]
         result[i, i, i] <- qux * (1 - qux) * (1 - 2 * qux)
     }
     return(result)
 }

 theta <- matrix(rnorm(d * 25), ncol = d)

 itself <- apply(theta, 1, cumfun)

 grad <- t(apply(theta, 1, moofun))
 all.equal(rowSums(grad), rep(1, nrow(grad)))

 fred <- function(theta) {
     stopifnot(length(theta) == d)
     result <- double(d)
     for (i in 1:d) {
         theta.plus <- theta
         theta.plus[i] <- theta.plus[i] + epsilon
         theta.minus <- theta
         theta.minus[i] <- theta.minus[i] - epsilon
         result[i] <- (cumfun(theta.plus) - cumfun(theta.minus)) / (2 * epsilon)
     }
     return(result)
 }
 fdgrad <- t(apply(theta, 1, fred))
 all.equal(grad, fdgrad)

 hess <- apply(theta, 1, voofun)
 hess <- array(as.vector(hess), dim = c(d, d, ncol(hess)))
 dim(hess)

 fred <- function(theta) {
     stopifnot(length(theta) == d)
     result <- matrix(NA, d, d)
     for (i in 1:d) {
         theta.plus <- theta
         theta.plus[i] <- theta.plus[i] + epsilon
         theta.minus <- theta
         theta.minus[i] <- theta.minus[i] - epsilon
         result[i, ] <- (moofun(theta.plus) - moofun(theta.minus)) /
             (2 * epsilon)
     }
     return(result)
 }

 fdhess <- apply(theta, 1, fred)
 fdhess <- array(as.vector(fdhess), dim = c(d, d, ncol(fdhess)))
 all.equal(hess, fdhess)

 thss <- apply(theta, 1, thoofun)
 thss <- array(as.vector(thss), dim = c(d, d, d, ncol(thss)))
 dim(thss)

 fred <- function(theta) {
     stopifnot(length(theta) == d)
     result <- array(NA, dim = c(d, d, d))
     for (i in 1:d) {
         theta.plus <- theta
         theta.plus[i] <- theta.plus[i] + epsilon
         theta.minus <- theta
         theta.minus[i] <- theta.minus[i] - epsilon
         result[i, , ] <- (voofun(theta.plus) - voofun(theta.minus)) /
             (2 * epsilon)
     }
     return(result)
 }
 fdthss <- apply(theta, 1, fred)
 fdthss <- array(as.vector(fdthss), dim = c(d, d, d, ncol(fdthss)))
 all.equal(thss, fdthss)

 cumfun <- function(theta) cumulant(theta, fam.multinomial(d))$zeroth
 moofun <- function(theta) cumulant(theta, fam.multinomial(d), deriv = 1)$first
 voofun <- function(theta) cumulant(theta, fam.multinomial(d), deriv = 2)$second
 thoofun <- function(theta) cumulant(theta, fam.multinomial(d), deriv = 3)$third

 all.equal(itself, apply(theta, 1, cumfun))
 all.equal(grad, t(apply(theta, 1, moofun)))
 chess <- apply(theta, 1, voofun)
 chess <- array(as.vector(chess), dim = c(d, d, ncol(chess)))
 all.equal(hess, chess)
 cthss <- apply(theta, 1, thoofun)
 cthss <- array(as.vector(cthss), dim = c(d, d, d, ncol(cthss)))
 all.equal(cthss, thss)

 foo <- apply(grad, 1, function(xi) link(xi, fam.multinomial(d),
     deriv = 0)$zeroth)
 bar <- apply(t(foo), 1, function(theta) cumulant(theta, fam.multinomial(d),
     deriv = 1)$first)
 all.equal(t(bar), grad)

 baz <- apply(grad, 1, function(xi) link(xi, fam.multinomial(d),
     deriv = 1)$first)
 baz <- array(as.vector(baz), dim = c(d, d, nrow(theta)))
 fdjackfun <- function(xi) {
     result <- matrix(NA, d, d)
     for (i in 1:d) {
         xi.plus <- xi
         xi.plus[i] <- xi[i] + epsilon
         xi.minus <- xi
         xi.minus[i] <- xi[i] - epsilon
         result[ , i] <- (link(xi.plus, fam.multinomial(d))$zeroth -
             link(xi.minus, fam.multinomial(d))$zeroth) / (2 * epsilon)
     }
     return(result)
 }
 qux <- apply(grad, 1, fdjackfun)
 qux <- array(as.vector(qux), dim = c(d, d, nrow(theta)))
 all.equal(baz, qux)

 ##### Normal location-scale #####

 d <- 2

 cumfun <- function(theta) {
     stopifnot(length(theta) == d)
     stopifnot(theta[2] < 0)
     - theta[1]^2 / (4 * theta[2]) + (1 / 2) * log(- 1 / (2 * theta[2]))
 }

 moofun <- function(theta) {
     stopifnot(length(theta) == d)
     stopifnot(theta[2] < 0)
     r1 <- (- theta[1] / (2 * theta[2]))
     r2 <- theta[1]^2 / (4 * theta[2]^2) - 1 / (2 * theta[2])
     return(c(r1, r2))
 }

 voofun <- function(theta) {
     stopifnot(length(theta) == d)
     stopifnot(theta[2] < 0)
     r11 <- (- 1 / (2 * theta[2]))
     r12 <- theta[1] / (2 * theta[2]^2)
     r22 <- (- theta[1]^2 / (2 * theta[2]^3)) + 1 / (2 * theta[2]^2)
     return(matrix(c(r11, r12, r12, r22), 2, 2))
 }

 thoofun <- function(theta) {
     stopifnot(length(theta) == d)
     stopifnot(theta[2] < 0)
     r111 <- 0
     r112 <- 1 / (2 * theta[2]^2)
     r122 <- (- theta[1] / theta[2]^3)
     r222 <- 3 * theta[1]^2 / (2 * theta[2]^4) - 1 / theta[2]^3
     return(array(c(r111, r112, r112, r122, r112, r122, r122, r222), dim = c(2, 2, 2)))
 }

 theta <- cbind(rnorm(25), - rexp(25))

 itself <- apply(theta, 1, cumfun)

 grad <- t(apply(theta, 1, moofun))

 fred <- function(theta) {
     stopifnot(length(theta) == d)
     result <- double(d)
     for (i in 1:d) {
         theta.plus <- theta
         theta.plus[i] <- theta.plus[i] + epsilon
         theta.minus <- theta
         theta.minus[i] <- theta.minus[i] - epsilon
         result[i] <- (cumfun(theta.plus) - cumfun(theta.minus)) / (2 * epsilon)
     }
     return(result)
 }
 fdgrad <- t(apply(theta, 1, fred))
 all.equal(grad, fdgrad)

 hess <- apply(theta, 1, voofun)
 hess <- array(as.vector(hess), dim = c(d, d, ncol(hess)))
 dim(hess)

 fred <- function(theta) {
     stopifnot(length(theta) == d)
     result <- matrix(NA, d, d)
     for (i in 1:d) {
         theta.plus <- theta
         theta.plus[i] <- theta.plus[i] + epsilon
         theta.minus <- theta
         theta.minus[i] <- theta.minus[i] - epsilon
         result[i, ] <- (moofun(theta.plus) - moofun(theta.minus)) /
             (2 * epsilon)
     }
     return(result)
 }

 fdhess <- apply(theta, 1, fred)
 fdhess <- array(as.vector(fdhess), dim = c(d, d, ncol(fdhess)))
 all.equal(hess, fdhess)

 thss <- apply(theta, 1, thoofun)
 thss <- array(as.vector(thss), dim = c(d, d, d, ncol(thss)))
 dim(thss)

 fred <- function(theta) {
     stopifnot(length(theta) == d)
     result <- array(NA, dim = c(d, d, d))
     for (i in 1:d) {
         theta.plus <- theta
         theta.plus[i] <- theta.plus[i] + epsilon
         theta.minus <- theta
         theta.minus[i] <- theta.minus[i] - epsilon
         result[i, , ] <- (voofun(theta.plus) - voofun(theta.minus)) /
             (2 * epsilon)
     }
     return(result)
 }

 fdthss <- apply(theta, 1, fred)
 fdthss <- array(as.vector(fdthss), dim = c(d, d, d, ncol(fdthss)))
 all.equal(thss, fdthss)

 cumfun <- function(theta) cumulant(theta, fam.normal.location.scale())$zeroth
 moofun <- function(theta) cumulant(theta, fam.normal.location.scale(),
     deriv = 1)$first
 voofun <- function(theta) cumulant(theta, fam.normal.location.scale(),
     deriv = 2)$second
 thoofun <- function(theta) cumulant(theta, fam.normal.location.scale(),
     deriv = 3)$third

 all.equal(itself, apply(theta, 1, cumfun))
 all.equal(grad, t(apply(theta, 1, moofun)))
 chess <- apply(theta, 1, voofun)
 chess <- array(as.vector(chess), dim = c(d, d, ncol(chess)))
 all.equal(hess, chess)
 cthss <- apply(theta, 1, thoofun)
 cthss <- array(as.vector(cthss), dim = c(d, d, d, ncol(cthss)))
 all.equal(cthss, thss)

 foo <- apply(grad, 1, function(xi) link(xi, fam.normal.location.scale(),
     deriv = 0)$zeroth)
 foo <- t(foo)
 all.equal(theta, foo)

 foo <- apply(grad, 1, function(xi) link(xi, fam.normal.location.scale(),
     deriv = 1)$first)
 foo <- array(as.vector(foo), dim = c(2, 2, nrow(theta)))
 is.ok <- as.logical(rep(0, nrow(theta)))
 for (i in 1:nrow(theta))
     is.ok[i] <- all.equal(foo[ , , i] %*% hess[ , , i], diag(2))
 all(is.ok)

