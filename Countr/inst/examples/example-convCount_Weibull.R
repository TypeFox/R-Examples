## user pwei
pwei_user <- function(tt, distP) {
    alpha <- exp(-log(distP[["scale"]]) / distP[["shape"]])
    pweibull(q = tt, scale = alpha, shape = distP[["shape"]],
             lower.tail = FALSE)
}


test_that("test all-probs convolution weibull count", {
    print("~~~~~~ Testing all-probs convolution weibull count ~~~~~~~~")
    x <- 0:10
    lambda <- 2.56
    p0 <- dpois(x, lambda)
    ll <- sum(dpois(x, lambda, TRUE))

    err <- 1e-6
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ probability ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ============================== allProbs =================================
    print("***** all-probs conv approach ...")
    distPars <- list(scale = lambda, shape = 1)
    pmat_bi <- dCount_conv_bi(x, distPars, "weibull", "direct",
                              nsteps = 200)
    pmat_user <- dCount_conv_user(x, distPars, c(1, 2), pwei_user, "direct",
                                  nsteps = 200)
    expect_less_than(max((pmat_bi- p0)^2 / p0), err)
    expect_less_than(max((pmat_user- p0)^2 / p0), err)

    ## ============================== naive =================================
    print("***** naive conv approach ...")
    pmat_bi <- dCount_conv_bi(x, distPars, "weibull", "naive",
                              nsteps = 200)
    pmat_user <- dCount_conv_user(x, distPars, c(1, 2), pwei_user, "naive",
                                  nsteps = 200)

    expect_less_than(max((pmat_bi- p0)^2 / p0), err)
    expect_less_than(max((pmat_user- p0)^2 / p0), err)

    ## ============================== dePril =================================
    print("***** dePril conv approach ...")
    pmat_bi <- dCount_conv_bi(x, distPars, "weibull", "dePril",
                              nsteps = 200)
    pmat_user <- dCount_conv_user(x, distPars, c(1, 2), pwei_user, "dePril",
                                  nsteps = 200)

    expect_less_than(max((pmat_bi- p0)^2 / p0), err)
    expect_less_than(max((pmat_user- p0)^2 / p0), err)

    ## ~~~~~~~~~~~~~~~~~~~ log-likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    distParsList <- lapply(seq(along = x), function(ind) distPars)
    extrapolParsList <- lapply(seq(along = x), function(ind) c(2, 1))

    ## ============================== allProbs =================================
    print("log-likehood allProbs Poisson ...")
    expect_equal(dCount_conv_loglik_bi(x, distParsList,
                                       "weibull", "direct", nsteps = 400),
                 ll, tolerance = 1e-4)
    expect_equal(dCount_conv_loglik_user(x, distParsList, extrapolParsList,
                                             pwei_user, "direct", nsteps = 400),
                 ll, tolerance = 1e-4)
    ## ============================== naive =================================
    print("log-likehood naive Poisson ...")
    expect_equal(dCount_conv_loglik_bi(x, distParsList,
                                       "weibull", "naive", nsteps = 400),
                 ll, tolerance = 1e-4)
    expect_equal(dCount_conv_loglik_user(x, distParsList, extrapolParsList,
                                             pwei_user, "naive", nsteps = 400),
                 ll, tolerance = 1e-4)
    ## ============================== dePril =================================
    print("log-likehood dePril Poisson ...")
    expect_equal(dCount_conv_loglik_bi(x, distParsList,
                                           "weibull", "dePril", nsteps = 400),
                 ll, tolerance = 1e-4)
    expect_equal(dCount_conv_loglik_user(x, distParsList, extrapolParsList,
                                             pwei_user, "dePril", nsteps = 400),
                 ll, tolerance = 1e-4)
})

test_that("test all-probs convolution weibull count -- McShane results", {
    print("~~~~~~ Testing all-probs convolution weibull count ---McShane results ~~~~~~~~")
                # res <- readRDS("McShane_paperResults.RDS")
    fn <- system.file("extdata", "McShane_paperResults.RDS", package = "Countr")
    res <- readRDS(fn)

    y <- res$y
    ## =========================== all-Probs =====================================
    print("***** all-probs approach ...")
    reswei <- res$weibullCountRes
    distPars <- list(scale = reswei$scale, shape = reswei$shape)
    distParsList <- lapply(seq(along = y), function(ind) distPars)
    extrapolParsList <- lapply(seq(along = y), function(ind) c(distPars$shape, 2))
    ll_bi <- dCount_conv_loglik_bi(y, distParsList,
                                       "weibull", "direct", nsteps = 200)
    ll_user <- dCount_conv_loglik_user(y, distParsList, extrapolParsList,
                                           pwei_user, "direct", nsteps = 200)

    expect_less_than(abs(ll_bi - reswei$loglik), 0.4)
    expect_less_than(abs(ll_user - reswei$loglik), 0.4)

    ## =========================== naive =====================================
    print("***** naive approach ...")
    reswei <- res$weibullCountRes
    distPars <- list(scale = reswei$scale, shape = reswei$shape)
    distParsList <- lapply(seq(along = y), function(ind) distPars)
    extrapolParsList <- lapply(seq(along = y), function(ind) c(distPars$shape, 2))
    ll_bi <- dCount_conv_loglik_bi(y, distParsList,
                                       "weibull", "naive", nsteps = 200)
    ll_user <- dCount_conv_loglik_user(y, distParsList, extrapolParsList,
                                           pwei_user, "naive", nsteps = 200)

    expect_less_than(abs(ll_bi - reswei$loglik), 0.4)
    expect_less_than(abs(ll_user - reswei$loglik), 0.4)

    ## =========================== dePril =====================================
    print("***** dePril approach ...")
    reswei <- res$weibullCountRes
    distPars <- list(scale = reswei$scale, shape = reswei$shape)
    distParsList <- lapply(seq(along = y), function(ind) distPars)
    extrapolParsList <- lapply(seq(along = y), function(ind) c(distPars$shape, 2))
    ll_bi <- dCount_conv_loglik_bi(y, distParsList,
                                   "weibull", "dePril", nsteps = 200)
    ll_user <- dCount_conv_loglik_user(y, distParsList, extrapolParsList,
                                       pwei_user, "dePril", nsteps = 200)

    expect_less_than(abs(ll_bi - reswei$loglik), 0.4)
    expect_less_than(abs(ll_user - reswei$loglik), 0.4)
})
