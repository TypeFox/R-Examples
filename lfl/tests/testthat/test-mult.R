test_that('mult', {
    set.seed(4523)
    x <- matrix(runif(24, -100, 100), ncol=6)
    y <- matrix(runif(18, -100, 100), nrow=6)
    res <- mult(x, y, function(xx, yy) sum(xx * yy))
    expect_true(is.numeric(res))
    expect_true(is.matrix(res))
    expect_equal(res, x %*% y)
})

