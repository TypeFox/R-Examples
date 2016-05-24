
 library(rcdd)

 set.seed(42)

 n <- 50
 d <- 5
 x <- matrix(rnorm(n * d), nrow = n)
 x <- d2q(x)

 foo <- cbind("0", cbind("1", x))
 # in order to get the same results as before save and restore .Random.seed
 # because sccd now calls R RNG
 save.Random.seed <- .Random.seed
 out <- scdd(foo, representation = "V")
 .Random.seed <- save.Random.seed

 l <- out$output[ , 1]
 b <- out$output[ , 2]
 a <- out$output[ , - c(1, 2)]
 a <- qxq(a, rep(-1, length(a)))

 unique(l)

 axb <- qmatmult(a, t(x))
 axb <- sweep(axb, 1, b, FUN = qmq)
 fred <- apply(axb, 2, function(foo) max(qsign(foo)))

 all(fred <= 0)
 ### points in interior
 sum(fred < 0)
 ### points on boundary
 sum(fred == 0)

 ### try on some new points

 y <- matrix(rnorm(2 * n * d), nrow = 2 * n)
 y <- d2q(y)

 ### REVISED DOWN TO HERE

 ayb <- qmatmult(a, t(y))
 ayb <- sweep(ayb, 1, b, FUN = qmq)
 fred <- apply(ayb, 2, function(foo) max(qsign(foo)))

 ### points in interior
 sum(fred < 0)
 ### points on boundary
 sum(fred == 0)
 ### points in exterior
 sum(fred > 0)

