test_that("test built-in distributions", {
    print("~~~~~~ Testing built-in distributions ~~~~~~~~")
    tt <- 2.5
    ## ============================== weibull =================================
    print("***** weibull ...")
    distP <- list(scale = 1.2, shape = 1.16)
    alpha <- exp(-log(distP[["scale"]]) / distP[["shape"]])
    expect_equal(pweibull(q = tt, scale = alpha, shape = distP[["shape"]],
                          lower.tail = FALSE),
                 surv(tt, distP, "weibull"),
                 tolerance = 1e-3)
    ## =============================== gamma ===================================
    print("***** gamma ...")
    distP <- list(shape = 0.5, rate = 1.0 / 0.7)
    expect_equal(pgamma(q = tt, rate = distP[["rate"]], shape = distP[["shape"]],
                        lower.tail = FALSE),
                 surv(tt, distP, "gamma"),
                 tolerance = 1e-3)
    ## =============================== gengamma ================================
    print("***** generalized gamma ...")
    distP <- list(mu = 0.5, sigma = 0.7, Q = 0.7)
    expect_equal(flexsurv::pgengamma(q = tt, mu = distP[["mu"]],
                                     sigma = distP[["sigma"]],
                                     Q = distP[["Q"]],
                                     lower.tail = FALSE),
                 surv(tt, distP, "gengamma"),
                 tolerance = 1e-3)
})
