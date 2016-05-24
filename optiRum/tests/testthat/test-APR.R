context("APR")

test_that("APR correctly produces values, no FV", {
    check <- 0.1766  # Taken from excel
    expect_true(round(APR(12, -10, 110), 4) == check)
    check <- c(0.1766, 0.0884)  # Taken from excel
    df <- data.frame(nper = c(12, 24), pmt = c(-10, -10), pv = c(110, 220))
    expect_true(identical(round(APR(df$nper, df$pmt, df$pv), 4), check))
})

test_that("APR correctly produces values, FV", {
    check <- 0.0895  # Taken from excel
    expect_true(round(APR(12, -10, 110, 5), 4) == check)
    check <- c(0.0895, 0.0674)  # Taken from excel
    df <- data.frame(nper = c(12, 24), pmt = c(-10, -10), pv = c(110, 220), fv = c(5, 5))
    expect_true(identical(round(APR(df$nper, df$pmt, df$pv, df$fv), 4), check))
})


test_that("APR errors given incorrect inputs", {
    expect_error(APR(0, -500, 3000))
    expect_error(APR(1, 500, 3000))
    expect_error(APR(1, -500, -3000))
    expect_error(APR("0", -500, 3000))
    expect_error(APR(1, "500", 3000))
    expect_error(APR(1, -500, "-3000"))
}) 
