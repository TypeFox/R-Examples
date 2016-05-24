source("helper-diversitree.R")

context("Unresolved clades (BiSSE + BiSSE-ness)")

make.branches.dtlik <- diversitree:::make.branches.dtlik
unresolv.b <- diversitree:::branches.unresolved.bisse
unresolv.n <- diversitree:::branches.unresolved.bisseness

control <- diversitree:::check.control.ode(list())
info.b <- diversitree:::make.info.bisse(NULL)
info.n <- diversitree:::make.info.bisseness(NULL)

branches.b <- make.branches.dtlik(info.b, control)
branches.n <- make.branches.dtlik(info.n, control)

tt <- seq(0, 5, length.out=21)
y0 <- c(0, 0, 1, 0)
y1 <- c(0, 0, 0, 1)
u0 <- c(as.list(data.frame(Nc=1L, k=0L, nsc=1L, len=tt)), list(nt.extra=20))
u1 <- c(as.list(data.frame(Nc=1L, k=1L, nsc=1L, len=tt)), list(nt.extra=20))

set.seed(1)
pars.b <- sort(runif(6))
pars.n <- c(sort(runif(6)), runif(4))

out.0.b.b <- branches.b(y0, tt, pars.b, 0, NA)
out.0.b.u <- unresolv.b(pars.b, u0)
expect_that(out.0.b.b, equals(out.0.b.u[2:3], check.attributes=FALSE))

out.1.b.b <- branches.b(y1, tt, pars.b, 0, NA)
out.1.b.u <- unresolv.b(pars.b, u1)
expect_that(out.1.b.b, equals(out.1.b.u[2:3], check.attributes=FALSE))

out.0.n.b <- branches.n(y0, tt, pars.n, 0, NA)
out.0.n.u <- unresolv.n(pars.n, u0)
expect_that(out.0.b.b, equals(out.0.b.u[2:3], check.attributes=FALSE))

out.1.n.b <- branches.n(y1, tt, pars.n, 0, NA)
out.1.n.u <- unresolv.n(pars.n, u1)
expect_that(out.1.b.b, equals(out.1.b.u[2:3], check.attributes=FALSE))
