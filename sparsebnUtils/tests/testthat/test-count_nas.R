context("count_nas")

test_that("count_nas behaves as expected", {
    expect_error(count_nas(rep(NA, 3)), "Input must be")

    ### Test a column of all NAs
    df.nrows <- 10
    df <- data.frame(x = rep(NA, df.nrows), y = runif(df.nrows))
    expect_equal(count_nas(df), nrow(df))

    ### Test no missing values
    df <- data.frame(x = rnorm(df.nrows), y = runif(df.nrows))
    expect_equal(count_nas(df), 0)

    ### Add some missing values at the bottom
    df[df.nrows + 1, ] <- c(NA, 1)
    df[df.nrows + 2, ] <- c(NA, NA)
    expect_equal(count_nas(df), 3)
})
