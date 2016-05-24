library(kappalab)

mu <- card.set.func(rnorm(9))
stopifnot(is.cardinal(mu) == is.cardinal(as.set.func(mu)))



