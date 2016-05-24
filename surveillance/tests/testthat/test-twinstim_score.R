context("Likelihood and score function of twinstim()")
## Note: derivatives of interaction functions are tested in separate files
##       we thus use relatively fast step functions here

data("imdepi")
model <- twinstim(
    endemic = addSeason2formula(~offset(log(popdensity)),
        S = 1, period = 365, timevar = "start"),
    epidemic = ~type,
    siaf = siaf.step(c(5, 20), maxRange = 50),
    tiaf = tiaf.step(2),
    data = imdepi,
    optim.args = NULL, verbose = FALSE
)
theta <- c("h.(Intercept)" = -20,
           "h.sin(2 * pi * start/365)" = 0.2, "h.cos(2 * pi * start/365)" = 0.3,
           "e.(Intercept)" = -10, "e.typeC" = -0.9,
           "e.siaf.1" = -1, "e.siaf.2" = -3, "e.tiaf.1" = -1)

test_that("likelihood is still the same", {
    expect_that(model$ll(theta), equals(-9610.68695991737))
})

test_that("score vector agrees with numerical approximation", {
    numsc <- if (surveillance.options("allExamples") && requireNamespace("numDeriv")) {
        numDeriv::grad(func = model$ll, x = theta)
    } else { # for faster --as-cran tests
        c(-365.19927878021, -29.3546236207476, -45.8139085706014,
          -88.5862997849202, -24.0808271983838, -54.0273836522059,
          -28.0233414216383, -74.5539641345285)
    }
    expect_that(model$sc(theta), equals(numsc))
})

## Note: twinstim() uses an estimate of the _expected_ Fisher information,
##       which does not necessarily agree with the negative Hessian of the ll
##       (it does asymptotically at the MLE)
## numfi <- -numDeriv::hessian(func = model$ll, x = theta)
## anafi <- model$fi(theta)


### now check with identity link for the epidemic predictor

model2 <- update(model, siaf = NULL, tiaf = NULL, epidemic = ~1, epilink = "log")
model2i <- update(model2, epilink = "identity")
theta2 <- theta2i <- theta[1:4]
theta2i["e.(Intercept)"] <- exp(theta2["e.(Intercept)"])

test_that("likelihoods with log-link and identity link are the same", {
    expect_that(model2$ll(theta2), equals(model2i$ll(theta2i)))
})

test_that("identity link score vector agrees with numerical approximation", {
    numsc <- if (surveillance.options("allExamples") && requireNamespace("numDeriv")) {
        numDeriv::grad(func = model2i$ll, x = theta2i)
    } else { # for faster --as-cran tests
        c(-679.706275919901, -91.0659401491325, -114.082117122738,
          -1532144485.45524)
    }
    expect_that(model2i$sc(theta2i), equals(numsc))
})
