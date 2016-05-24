library(kappalab)

a <- Mobius.game(c(0,runif(15)),4,4)
mu <- zeta(a)
 
f <- c(0.2,0.3,0.1,0.7)
stopifnot(abs(Choquet.integral(a,f) - Choquet.integral(mu,f)) < 1e-6)
stopifnot(abs(Sugeno.integral(a,f) - Sugeno.integral(mu,f)) < 1e-6)
f <- c(0.2,-0.3,0.1,-0.7)
stopifnot(abs(Sipos.integral(a,f) - Sipos.integral(mu,f)) < 1e-6)



