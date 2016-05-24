library(kappalab)

mu <- capacity(c(0:13,13,13)/13)
a <- Mobius(mu)
stopifnot(abs(normalize(normalize(mu))@data - normalize(mu)@data) < 1e-6)
stopifnot(abs(mu@data - conjugate(conjugate(mu))@data) < 1e-6)
stopifnot(abs(mu@data - zeta(Mobius(mu))@data) < 1e-6)
stopifnot(abs(a@data - Mobius(zeta(a))@data) < 1e-6)
stopifnot(abs(veto(mu) - veto(a)) < 1e-6)
stopifnot(abs(favor(mu) - favor(a)) < 1e-6)
stopifnot(abs(orness(mu) - orness(a)) < 1e-6)
stopifnot(abs(variance(mu) - variance(a)) < 1e-6)
stopifnot(abs(entropy(mu) - entropy(a)) < 1e-6)



