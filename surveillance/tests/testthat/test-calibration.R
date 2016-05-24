context("Calibration tests for Poisson or NegBin predictions")

mu <- c(0.1, 1, 3, 6, pi, 100)
size1 <- 0.5
size2 <- c(0.1, 0.1, 10, 10, 100, 100)
##set.seed(2); y <- rnbinom(length(mu), mu = mu, size = size1)
y <- c(0, 0, 2, 14, 5, 63)

zExpected <- rbind(
    dss = c(P = 6.07760977730636, NB1 = -0.468561113465647, NB2 = 2.81071829075294),
    logs = c(P = 5.95533908528874, NB1 = 0.403872251419915, NB2 = 2.77090543018323),
    rps = c(P = 4.45647234878906, NB1 = -0.437254253267393, NB2 = 2.57223607389215)
    )

delta <- 1e-4 #sqrt(.Machine$double.eps)

for (score in rownames(zExpected)) {
    .zExpected <- zExpected[score, , drop = TRUE]
    ## if package "gsl" is not available, rps_EV is less accurate
    tol_equal <- if (score == "rps" && !requireNamespace("gsl", quietly = TRUE))
                     1e-4 else .Machine$double.eps^0.5
    test_that(paste0("still the same z-statistics with ", score), {
        ## Poisson predictions
        zP <- calibrationTest(y, mu, which = score, tolerance = delta)$statistic
        expect_equal(zP, .zExpected["P"], check.attributes = FALSE,
                     tolerance = tol_equal)
        ## NegBin predictions with common size parameter
        zNB1 <- calibrationTest(y, mu, size1, which = score, tolerance = delta)$statistic
        expect_equal(zNB1, .zExpected["NB1"], check.attributes = FALSE,
                     tolerance = tol_equal)
        ## NegBin predictions with varying size parameter
        zNB2 <- calibrationTest(y, mu, size2, which = score, tolerance = delta)$statistic
        expect_equal(zNB2, .zExpected["NB2"], check.attributes = FALSE,
                     tolerance = tol_equal)
    })
}
