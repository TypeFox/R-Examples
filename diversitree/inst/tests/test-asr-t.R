source("helper-diversitree.R")

control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)

context("ASR (time-varing BiSSE)")

## Start with a simple tree evolved under a BiSSE with all rates
## asymmetric:
pars <- c(.1, .2, .03, .06, .01, .02)
set.seed(3)
phy <- trees(pars, "bisse", max.taxa=50, max.t=Inf, x0=0)[[1]]

## BiSSE ancestral state reconstructions under the ML model
lik.b <- make.bisse(phy, phy$tip.state)
lik.m <- make.musse(phy, phy$tip.state+1, 2)

## fit <- find.mle(lik, pars, method="subplex")
## p.b.ml <- dput(coef(fit))
p.b <- pars
p.b.ml <- structure(c(0.0840286576638916, 0.211042245168045,
                      0.0607003387416787, 3.34756696282173e-06,
                      0.00412941736207307, 0.0236232224532126),
                    .Names=c("lambda0", "lambda1", "mu0", "mu1", "q01",
                      "q10"))

st.b <- asr.marginal(lik.b, p.b)
st.m <- asr.marginal(lik.m, p.b)
expect_that(st.b, equals(st.m))

## Now, time dependence:
functions <- rep(c("linear.t", "constant.t"), c(2, 4))
lik.b.t.0 <- make.bisse.t(phy, phy$tip.state, functions)
lik.b.t.d <- make.bisse.t(phy, phy$tip.state, functions, control=control.d)
lik.b.t.g <- make.bisse.t(phy, phy$tip.state, functions, control=control.g)
lik.b.t.G <- make.bisse.t(phy, phy$tip.state, functions, control=control.G)

lik.m.t.0 <- make.musse.t(phy, phy$tip.state+1, 2, functions)
lik.m.t.d <- make.musse.t(phy, phy$tip.state+1, 2, functions, control=control.d)
lik.m.t.g <- make.musse.t(phy, phy$tip.state+1, 2, functions, control=control.g)
lik.m.t.G <- make.musse.t(phy, phy$tip.state+1, 2, functions, control=control.G)

## Time-independent parameters for time models
p.b.t <- c(rbind(p.b[1:2], 0), p.b[-(1:2)])
names(p.b.t) <- argnames(lik.b.t.0)

## Time-dependent parameters for time models
p.b.tt <- p.b.t
p.b.tt[sprintf("lambda%d.m", 0:1)] <- .002

## Sanity check -- make sure that the likelihoods agree.
## ll.t <- lik.b.t.0(p.b.t)
ll.t <- -157.302523694275
expect_that(lik.b.t.0(p.b.t), equals(ll.t))
expect_that(lik.b.t.d(p.b.t), equals(ll.t, tolerance=1e-7))
expect_that(lik.b.t.g(p.b.t), equals(ll.t))
expect_that(lik.b.t.G(p.b.t), equals(ll.t))
expect_that(lik.m.t.0(p.b.t), equals(ll.t))
expect_that(lik.m.t.d(p.b.t), equals(ll.t, tolerance=1e-7))
expect_that(lik.m.t.g(p.b.t), equals(ll.t))
expect_that(lik.m.t.G(p.b.t), equals(ll.t))

## And that they agree under a time-varying model
## ll.tt <- lik.b.t.0(p.b.tt)
ll.tt <- -161.035371378992
expect_that(lik.b.t.0(p.b.tt), equals(ll.tt))
expect_that(lik.b.t.d(p.b.tt), equals(ll.tt, tolerance=1e-7))
expect_that(lik.b.t.g(p.b.tt), equals(ll.tt))
expect_that(lik.b.t.G(p.b.tt), equals(ll.tt))
expect_that(lik.m.t.0(p.b.tt), equals(ll.tt))
expect_that(lik.m.t.d(p.b.tt), equals(ll.tt, tolerance=1e-7))
expect_that(lik.m.t.g(p.b.tt), equals(ll.tt))
expect_that(lik.m.t.G(p.b.tt), equals(ll.tt))

## Check the ancestral states:
st.b.t.0 <- asr.marginal(lik.b.t.0, p.b.t)
st.b.t.d <- asr.marginal(lik.b.t.d, p.b.t)
st.b.t.g <- asr.marginal(lik.b.t.g, p.b.t)
st.b.t.G <- asr.marginal(lik.b.t.G, p.b.t)

st.m.t.0 <- asr.marginal(lik.m.t.0, p.b.t)
st.m.t.d <- asr.marginal(lik.m.t.d, p.b.t)
st.m.t.g <- asr.marginal(lik.m.t.g, p.b.t)
st.m.t.G <- asr.marginal(lik.m.t.G, p.b.t)

expect_that(st.b.t.0, is_identical_to(st.b))
expect_that(st.m.t.0, is_identical_to(st.m))

expect_that(st.b.t.d, equals(st.b.t.0, tolerance=1e-7))
expect_that(st.b.t.g, equals(st.b.t.0))
expect_that(st.b.t.G, equals(st.b.t.0))
expect_that(st.m.t.d, equals(st.m.t.0, tolerance=1e-7))
expect_that(st.m.t.g, equals(st.m.t.0))
expect_that(st.m.t.G, equals(st.m.t.0))

## And with actual time varying...
st.b.tt.0 <- asr.marginal(lik.b.t.0, p.b.tt)
st.b.tt.d <- asr.marginal(lik.b.t.d, p.b.tt)
st.b.tt.g <- asr.marginal(lik.b.t.g, p.b.tt)
st.b.tt.G <- asr.marginal(lik.b.t.G, p.b.tt)

st.m.tt.0 <- asr.marginal(lik.m.t.0, p.b.tt)
st.m.tt.d <- asr.marginal(lik.m.t.d, p.b.tt)
st.m.tt.g <- asr.marginal(lik.m.t.g, p.b.tt)
st.m.tt.G <- asr.marginal(lik.m.t.G, p.b.tt)

## Should differ from time invariant version
expect_that(isTRUE(all.equal(st.b.tt.0, st.b.t.0)), is_false())

## Some point measurements
expect_that(var(st.b.tt.0[1, ]), equals(0.144781638970694))
expect_that(st.b.tt.0[1, 3], equals(0.620067231484044))

## Internally consistent:
expect_that(st.b.t.d, equals(st.b.t.0, tolerance=1e-7))
expect_that(st.b.t.g, equals(st.b.t.0))
expect_that(st.b.t.G, equals(st.b.t.0))
expect_that(st.m.t.0, equals(st.m.t.0))
expect_that(st.m.t.d, equals(st.m.t.0, tolerance=1e-7))
expect_that(st.m.t.g, equals(st.m.t.0))
expect_that(st.m.t.G, equals(st.m.t.0))

context("ASR (time-varing MuSSE)")

pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32

set.seed(2)
phy <- tree.musse(pars, 30, x0=1)
states <- phy$tip.state

functions <- rep(c("linear.t", "constant.t"), c(3, 9))

lik.m <- make.musse(phy, states, 3)
lik.t.0 <- make.musse.t(phy, states, 3, functions)
lik.t.d <- make.musse.t(phy, states, 3, functions, control=control.d)
lik.t.g <- make.musse.t(phy, states, 3, functions, control=control.g)
lik.t.G <- make.musse.t(phy, states, 3, functions, control=control.G)

p <- starting.point.musse(phy, 3)
p.t <- c(rbind(p[1:3], 0), p[-(1:3)])
names(p.t) <- argnames(lik.t.0)

ll <- lik.m(p)
expect_that(lik.t.0(p.t), equals(ll))
expect_that(lik.t.d(p.t), equals(ll))
expect_that(lik.t.g(p.t), equals(ll))
expect_that(lik.t.G(p.t), equals(ll))

## Ancestral states:
st.m <- asr.marginal(lik.m, p)
st.t.0 <- asr.marginal(lik.t.0, p.t)
st.t.d <- asr.marginal(lik.t.d, p.t)
st.t.g <- asr.marginal(lik.t.g, p.t)
st.t.G <- asr.marginal(lik.t.G, p.t)

expect_that(st.t.0, is_identical_to(st.m))
expect_that(st.t.d, equals(st.m, tolerance=3e-7))
expect_that(st.t.g, equals(st.m))
expect_that(st.t.G, equals(st.m))

p.tt <- p.t
p.tt[sprintf("lambda%d.m", 1:3)] <- c(.02, .015, .01)

st.tt.0 <- asr.marginal(lik.t.0, p.tt)
st.tt.d <- asr.marginal(lik.t.d, p.tt)
st.tt.g <- asr.marginal(lik.t.g, p.tt)
st.tt.G <- asr.marginal(lik.t.G, p.tt)

expect_that(rowSums(st.tt.0),
            equals(c(11.4584448192075, 13.3218259039656, 
                     4.21972927682697)))
expect_that(st.tt.d, equals(st.tt.0, tolerance=2e-6))
expect_that(st.tt.g, equals(st.tt.0))
expect_that(st.tt.G, equals(st.tt.0))
