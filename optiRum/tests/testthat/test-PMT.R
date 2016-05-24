context("PMT")

test_that("PMT correctly produces values", {
    check <- -440.29  # Taken from excel
    expect_true(PMT(0.1, 12, 3000) == check)
    check <- c(-440.29, -202.55)  # Taken from excel
    df <- data.frame(rate = c(0.1, 0.2), nper = c(12, 24), pv = c(3000, 1000))
    expect_true(identical(PMT(df$rate, df$nper, df$pv), check))
})

test_that("PMT errors given incorrect inputs", {
    expect_error(PMT(0, 12, 3000))
    expect_error(PMT(0.1, 0, 3000))
    expect_error(PMT(0.1, 12, -3000))
    expect_error(PMT("0", 12, 3000))
    expect_error(PMT(0.1, "0", 3000))
    expect_error(PMT(0.1, 12, "-3000"))
}) 
