pwei_user <- function(tt, distP) {
    alpha <- exp(-log(distP[["scale"]]) / distP[["shape"]])
    pweibull(q = tt, scale = alpha, shape = distP[["shape"]],
             lower.tail = FALSE)
}


test_that("test ev convolution weibull count", {
    print("~~~~~~ Testing ev convolution Poisson count ~~~~~~~~")
    lambda <- 2.56
    beta <- 1

    distPars <- list(scale = lambda, shape = beta)
    evbi <- evCount_conv_bi(20, distPars, dist = "weibull")
    evu <- evCount_conv_user(20, distPars, c(2, 2), pwei_user, "dePril")

    expect_equal(evbi[["ExpectedValue"]], lambda, tolerance = 0.05)
    expect_equal(evu[["ExpectedValue"]], lambda, tolerance = 0.05)
    expect_equal(evbi[["Variance"]], lambda, tolerance = 0.05)
    expect_equal(evu[["Variance"]], lambda, tolerance = 0.05)

    print("~~~~~~ Testing ev convolution weibull count ~~~~~~~~")
    lambda <- 2.56
    beta <- 1.35

    distPars <- list(scale = lambda, shape = beta)
    evbi <- evCount_conv_bi(20, distPars, dist = "weibull")
    evu <- evCount_conv_user(20, distPars, c(2.35, 2), pwei_user, "dePril")

    x <- 1:20
    px <- dCount_conv_bi(x, distPars, "weibull", "dePril",
                         nsteps = 100)
    ev <- sum(x * px)
    var <- sum(x^2 * px) - ev^2
    
    expect_equal(evbi[["ExpectedValue"]], ev, tolerance = 0.05)
    expect_equal(evu[["ExpectedValue"]], ev, tolerance = 0.05)
    expect_equal(evbi[["Variance"]], var, tolerance = 0.05)
    expect_equal(evu[["Variance"]], var, tolerance = 0.05)

})
