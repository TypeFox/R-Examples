library(kappalab)

mu <- game(c(0,runif(31)))
stopifnot(abs(mu@data - conjugate(conjugate(mu))@data) < 1e-6)

mu <- capacity(c(0,rep(1,15)))
stopifnot(abs(mu@data - conjugate(conjugate(mu))@data) < 1e-6)

mu <- card.game(c(0,rnorm(17)))
stopifnot(abs(mu@data - conjugate(conjugate(mu))@data) < 1e-6)




