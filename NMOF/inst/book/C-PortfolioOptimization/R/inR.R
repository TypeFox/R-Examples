# inR.R -- version 2010-12-24
ns <- 100L; na <- 10L
R <- array(rnorm(ns * na), dim = c(ns, ns))
R1 <- crossprod(R)/ns 
R2 <- ((ns - 1)/ns)*cov(R) + outer(colMeans(R), colMeans(R))
identical(R1,R2); stopifnot(all.equal(R1,R2))

