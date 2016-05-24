context("Temporal interaction functions")

test_that("Step kernel of a single type agrees with numerical approximations",
     {
         steptiaf <- tiaf.step(c(7,20), maxRange=25, nTypes=1)
         logvals <- log(c(1.2,0.2))
         ##curve(steptiaf$g(x, logvals), 0, 30, n=301)

         ## check G
         expect_equal(steptiaf$G(0:30, logvals), sapply(0:30, function (upper) {
             integrate(steptiaf$g, 0, upper, logvals, rel.tol=1e-8)$value
         }), tolerance=1e-8)

         ## check deriv
         if (requireNamespace("maxLik", quietly = TRUE)) {
             expect_that(maxLik::compareDerivatives(
                 f = function(pars, x) steptiaf$g(x, pars),
                 grad = function(pars, x) steptiaf$deriv(x, pars),
                 t0 = logvals, x = c(0.5,2,5,7,10,15,20,25,30),
                 print = FALSE)$maxRelDiffGrad,
                         is_less_than(1e-8))
         }

         ## check Deriv
         for (paridx in seq_along(logvals))
             expect_equal(steptiaf$Deriv(0:30, logvals)[,paridx],
                          sapply(0:30, function (upper) integrate(
                              function(...) steptiaf$deriv(...)[,paridx],
                              0, upper, logvals, rel.tol=1e-6)$value),
                          tolerance=1e-6,
                          label=paste0("steptiaf$Deriv()[,",paridx,"]"),
                          expected.label="integrate() approximation")
                                 
     })

test_that("Step kernel with maxRange=eps.t is equivalent to maxRange=Inf",
      {
          data("imdepi", package="surveillance")
          imdfit_steptiafInf <- twinstim(
              endemic  = ~offset(log(popdensity)) + I(start/365 - 3.5),
              epidemic = ~1, siaf = siaf.constant(),
              tiaf = tiaf.step(c(7,20), maxRange=Inf),
              data = imdepi, optim.args = NULL, verbose = FALSE)
          maxepst <- max(imdepi$events$eps.t)
          imdfit_steptiaf30 <- update(
              imdfit_steptiafInf,
              tiaf = tiaf.step(c(7,20), maxRange=maxepst+0.1))
          coefs <- c(-20, -0.05, -15, -0.5, 0.2, -1)
          stopifnot(identical(imdfit_steptiafInf$ll(coefs),
                              imdfit_steptiaf30$ll(coefs)),
                    identical(imdfit_steptiafInf$sc(coefs),
                              imdfit_steptiaf30$sc(coefs)))
      })
