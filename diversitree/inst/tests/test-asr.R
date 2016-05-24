source("helper-diversitree.R")

context("ASR (BiSSE)")

## Simple BiSSE ASR:
pars <- c(.1, .2, .03, .06, .01, .02)
set.seed(3)
phy <- trees(pars, "bisse", max.taxa=50, max.t=Inf, x0=0)[[1]]

lik <- make.bisse(phy, phy$tip.state)
p <- c(0.084, 0.211, 0.061, 0, 0.004, 0.024)
st <- asr.marginal(lik, p)

expect_that(var(st[1, ]), equals(0.145473202035642))
expect_that(st[1, 3], equals(0.586896762121772))

set.seed(1)
p2 <- p + runif(length(p), 0, .1)
st2 <- asr.marginal(lik, p2)
expect_that(var(st2[1, ]), equals(0.0870348628276599))
expect_that(st2[1, 3], equals(0.064462088171836))

lik.m <- make.musse(phy, phy$tip.state+1, 2)
expect_that(asr.marginal(lik.m, p2), equals(st2))

lik.m <- make.mk2(phy, phy$tip.state)
fit.m <- find.mle(lik.m, pars[5:6], method="subplex")
st.m <- asr.marginal(lik.m, coef(fit.m))

st.id <- asr.marginal(lik, c(.1, .1, .03, .03, coef(fit.m)))
expect_that(st.id, equals(st.m, tolerance=1e-7))

set.seed(1)
st.id2 <- asr.marginal(lik, c(rep(runif(2, 0, .1), each=2), coef(fit.m)))
expect_that(st.id2, equals(st.m, tolerance=1e-7))

control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)

lik.b.d <- make.bisse(phy, phy$tip.state, control=control.d)
lik.b.g <- make.bisse(phy, phy$tip.state, control=control.g)
lik.b.G <- make.bisse(phy, phy$tip.state, control=control.G)

## The non-compiled versions are quite slow.
expect_that(asr.marginal(lik.b.d, p), equals(st, tolerance=5e-7))
expect_that(asr.marginal(lik.b.g, p), equals(st))
expect_that(asr.marginal(lik.b.G, p), equals(st))

expect_that(asr.marginal(lik.b.d, p2), equals(st2, tolerance=5e-7))
expect_that(asr.marginal(lik.b.g, p2), equals(st2))
expect_that(asr.marginal(lik.b.G, p2), equals(st2))

context("ASR (MuSSE)")
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
set.seed(2)
phy <- tree.musse(pars, 30, x0=1)

states <- phy$tip.state  
lik <- make.musse(phy, states, 3)

set.seed(1)
p <- pars + runif(length(pars), 0, .2)

expect_that(lik(p), equals(-118.952881268597))
st <- asr.marginal(lik, p)

expect_that(var(st[1, ]), equals(0.0521682452753316))
expect_that(st[1, 3], equals(0.0518022277082024))

lik.m <- make.mkn(phy, states, 3)
p2 <- p
p2[1:3] <- .1
p2[4:6] <- .03
p2.m <- p2[-(1:6)]

expect_that(lik(p2), equals(-119.15072625137572))
expect_that(lik.m(p2.m), equals(-28.64702621988359))

st1 <- asr.marginal(lik, p2)
st2 <- asr.marginal(lik.m, p2.m)
expect_that(st1, equals(st2, tolerance=5e-7))
