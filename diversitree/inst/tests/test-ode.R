## Testing for the ODE bits.
source("helper-diversitree.R")

## All testing is done with the BiSSE equations, just because they are
## the equations I've looked at most closely over the years.  In
## theory, this would be fairly easy to change to use to any model,
## provided parameter lengths and suitable parameters could be
## specified.
context("Solving ODEs")

derivs <- diversitree:::derivs.bisse
derivs.for.deSolve <- diversitree:::derivs.for.deSolve

## Initial conditions corresponding to a tip in state 0:
y <- c(0, 0, 1, 0)

## Vector of parameters
pars <- c(.1, .2, .03, .04, .01, .02)

## Vector of times to report at
tt <- seq(0, 30, length.out=101)

check.control.ode <- diversitree:::check.control.ode
control.gslode.C <- check.control.ode(list(backend="gslode",
                                           compiled=TRUE))
control.gslode.R <- check.control.ode(list(backend="gslode",
                                           compiled=FALSE))
control.deSolve <- check.control.ode(list(backend="deSolve"))

## Run the calculations with deSolve's lsoda (no clever tricks)
res.ref <- deSolve::lsoda(y, tt, derivs.for.deSolve(derivs), pars,
                          atol=control.deSolve$tol, rtol=control.deSolve$tol)
## Convert this to the format that we expect (drop time column and
## transpose data).
res.ref <- unname(t(res.ref[-1,-1,drop=FALSE]))

## Now build an "info" list:
info <-
  list(name="bisse", np=6, ny=4, idx.d=3:4, derivs=derivs,
       argnames=c("la0", "la1", "mu0", "mu1", "q01", "q10"))

## And build the ODE
make.ode <- diversitree:::make.ode
ode.deSolve  <- make.ode(info, control.deSolve)
ode.gslode.R <- make.ode(info, control.gslode.R)
ode.gslode.C <- make.ode(info, control.gslode.C)

res.deSolve  <- ode.deSolve(y,  tt, pars)
res.gslode.R <- ode.gslode.R(y, tt, pars)
res.gslode.C <- ode.gslode.C(y, tt, pars)

## The deSolve versions (ref and gslode.R) should be identical, as
## should the GslOde versions.  However, there should be numerical
## differences between the approaches.
expect_that(res.deSolve, is_identical_to(res.ref))

expect_that(identical(res.gslode.R, res.ref), is_false())
expect_that(res.gslode.R, equals(res.ref))
expect_that(res.gslode.C, equals(res.ref))
expect_that(res.gslode.R, is_identical_to(res.gslode.C))

## Now, try a time dependent model.
functions <- rep(c("linear.t", "constant.t"), c(2, 4))
names(functions) <- info$argnames

## I should expose this, I think (both?)
make.time.machine <- diversitree:::make.time.machine
make.derivs.t <- diversitree:::make.derivs.t

## TODO: the time range basically serves no purpose now.
tm <- make.time.machine(functions, range(tt))

## We have to set the parameters of the time machine before using it,
## though this is taken care of in the cases (add test)
pars.t <- c(pars[1], 0, pars[2], 0, pars[3:6])
tm$set(pars.t)

## Here is a time-dependent BiSSE derivative function.  Note that it
## ignores the parameter vector, accepting changes only through
## changes to tm.
derivs.t <- make.derivs.t(derivs, tm)

res.t.ref <- deSolve::lsoda(y, tt, derivs.for.deSolve(derivs.t), pars,
                            atol=control.deSolve$tol, rtol=control.deSolve$tol)
res.t.ref <- unname(t(res.t.ref[-1,-1,drop=FALSE]))

## Quick check that this agrees with the original reference set, as
## these parameters don't cause any actual time variation.
expect_that(res.t.ref, is_identical_to(res.ref))

## Then build time-varying functions through the usual make.ode
## interface.
info.t <- diversitree:::update.info.t(info, functions, range(tt))

ode.t.deSolve   <- make.ode(info.t, control.deSolve)
ode.t.gslode.R  <- make.ode(info.t, control.gslode.R)
ode.t.gslode.C  <- make.ode(info.t, control.gslode.C)

## OK -- this should work, but seems to fail, the same way for all
## three cases.
res.t.deSolve   <- ode.t.deSolve(y,  tt, pars.t)
res.t.gslode.R  <- ode.t.gslode.R(y, tt, pars.t)
res.t.gslode.C  <- ode.t.gslode.C(y, tt, pars.t)

## As above, the deSolve versions (ref and gslode.R) should be
## identical, as should the GslOde versions.  However, there should be
## numerical differences between the approaches.
expect_that(res.t.deSolve, is_identical_to(res.t.ref))

expect_that(identical(res.t.gslode.R, res.t.ref), is_false())
expect_that(res.t.gslode.R, equals(res.t.ref))
expect_that(res.t.gslode.C, equals(res.t.ref))
expect_that(res.t.gslode.R, is_identical_to(res.t.gslode.C))

## Snsure that the time machine parameters really are used.
tm$set(rep(0, 8))
res.t.deSolve   <- ode.t.deSolve(y,  tt, pars.t)
tm$set(rep(0, 8))
res.t.gslode.R  <- ode.t.gslode.R(y, tt, pars.t)
tm$set(rep(0, 8))
res.t.gslode.C  <- ode.t.gslode.C(y, tt, pars.t)

expect_that(res.t.deSolve, is_identical_to(res.t.ref))
expect_that(identical(res.t.gslode.R, res.t.ref), is_false())
expect_that(res.t.gslode.R, equals(res.t.ref))
expect_that(res.t.gslode.C, equals(res.t.ref))
expect_that(res.t.gslode.R, is_identical_to(res.t.gslode.C))

## Quick check to ensure that the compiled code really is used, by
## demanding a 10x speedup.
t.R <- system.time(ode.t.gslode.R(y, tt, pars.t))[["elapsed"]]
expect_that(ode.t.gslode.C(y, tt, pars.t),
            takes_less_than(t.R / 10))

## Or, more mysteriously, we could do:
## inherits(environment(ode.t.gslode.C)$ode, "Rcpp_GslOdeTime")
## inherits(environment(ode.t.gslode.R)$ode, "Rcpp_GslOdeR")

## Now, change the parameters, so that there is an effect of time.
pars.t[c(2,4)] <- .01
tm$set(pars.t)

## Rerun the above code
res.t.ref <- deSolve::lsoda(y, tt, derivs.for.deSolve(derivs.t), pars,
                            atol=control.deSolve$tol, rtol=control.deSolve$tol)
res.t.ref <- unname(t(res.t.ref[-1,-1,drop=FALSE]))

## Check that the reference example differs:
expect_that(isTRUE(all.equal(res.t.ref, res.ref)), is_false())

## The parameter 'pars.t' is actually ignored by the two .R versions,
## though the gslode.R one is using the wrong length.  Something
## should be done there.
res.t.deSolve   <- ode.t.deSolve(y, tt, pars.t)
res.t.gslode.R  <- ode.t.gslode.R(y, tt, pars.t)
res.t.gslode.C  <- ode.t.gslode.C(y, tt, pars.t)

expect_that(res.t.deSolve, is_identical_to(res.t.ref))
expect_that(identical(res.t.gslode.R, res.t.ref), is_false())
expect_that(res.t.gslode.R, equals(res.t.ref))
expect_that(res.t.gslode.C, equals(res.t.ref))
expect_that(res.t.gslode.R, is_identical_to(res.t.gslode.C))
