test_that('rmse', {
    f <- 1:10 + 0.2
    r <- 1:10

    expect_equal(rmse(f, r), 0.2)
})
