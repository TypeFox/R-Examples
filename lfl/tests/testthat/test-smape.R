test_that('smape', {
    f <- 1:10 + 0.2
    r <- 1:10

    expect_equal(smape(f, r), mean(0.2 / ((f+r) / 2)))
})
