context("PV")

test_that("PV correctly produces values, no fv", {
    check <- 68.14  # Taken from excel
    expect_true(PV(0.1, 12, -10) == check)
    check <- c(68.14, 134.77)  # Taken from excel
    df <- data.frame(rate = c(0.1, 0.1), nper = c(12, 24), pmt = c(-10, -15))
    expect_true(identical(PV(df$rate, df$nper, df$pmt), check))
})

test_that("PV correctly produces values, fv", {
    check <- 66.54  # Taken from excel
    expect_true(PV(0.1, 12, -10, 5) == check)
    check <- c(66.54, 134.26)  # Taken from excel
    df <- data.frame(rate = c(0.1, 0.1), nper = c(12, 24), pmt = c(-10, -15), fv = c(5, 5))
    expect_true(identical(PV(df$rate, df$nper, df$pmt, df$fv), check))
})


test_that("PV errors given incorrect inputs", {
    expect_error(PV(0, 12, -10))
    expect_error(PV(0.1, 0, -10))
    expect_error(PV(0.1, 12, 10))
    expect_error(PV("0", 12, -10))
    expect_error(PV(0.1, "12", -10))
    expect_error(PV(0.1, 12, "-10"))
    expect_error(PV(0.1, 12, -10, "5"))
}) 
