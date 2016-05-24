test_that('mase', {
    f <- 1:10 + 0.2
    r <- 1:10

    expect_equal(mase(f, r), 0.2 * 10 / ((10/9) * 9))
})
