library(kappalab)

mu <- card.game(c(0,rnorm(6)))
mu2 <- as.game(mu)
a <- Mobius(mu2)
stopifnot(abs(mu@data - as.card.game(as.game(mu))@data) < 1e-6)

f <- c(0.1,0.9,0.3,0.8,0.1,0.9)
stopifnot(abs(Choquet.integral(mu,f) - Choquet.integral(mu2,f)) < 1e-6)
stopifnot(abs(Sugeno.integral(mu,f) - Sugeno.integral(mu2,f)) < 1e-6)
stopifnot(abs(Choquet.integral(a,f) - Choquet.integral(mu2,f)) < 1e-6)
stopifnot(abs(Sugeno.integral(a,f) - Sugeno.integral(mu2,f)) < 1e-6)
stopifnot(abs(Choquet.integral(a,f) - Choquet.integral(mu,f)) < 1e-6)
stopifnot(abs(Sugeno.integral(a,f) - Sugeno.integral(mu,f)) < 1e-6)

f <- c(0.1,-0.9,-0.3,0.8,0.1,0.9)
stopifnot(abs(Sipos.integral(mu,f) - Sipos.integral(mu2,f)) < 1e-6)
stopifnot(abs(Sipos.integral(a,f) - Sipos.integral(mu2,f)) < 1e-6)
stopifnot(abs(Sipos.integral(a,f) - Sipos.integral(mu,f)) < 1e-6)
