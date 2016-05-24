test_that('sel.fsets', {
    orig <- matrix(runif(100), ncol=10)
    x <- fsets(orig,
               vars=1:10,
               specs=matrix(0, nrow=10, ncol=10))

    res <- sel(x, 1:5)
    expect_equal(as.matrix(res), orig[1:5, ])
    expect_equal(vars(res), vars(x))
    expect_equal(specs(res), specs(x))

    res <- sel(x, -3)
    expect_equal(as.matrix(res), orig[-3, ])
    expect_equal(vars(res), vars(x))
    expect_equal(specs(res), specs(x))

    res <- sel(x, 3)
    expect_equal(as.matrix(res), orig[3, , drop=FALSE])
    expect_equal(vars(res), vars(x))
    expect_equal(specs(res), specs(x))

    res <- sel(x, , 2:4)
    expect_equal(as.matrix(res), orig[, 2:4])
    expect_equal(vars(res), vars(x)[2:4])
    expect_equal(specs(res), specs(x)[2:4, 2:4])

    res <- sel(x, 1:5, 2:4)
    expect_equal(as.matrix(res), orig[1:5, 2:4])
    expect_equal(vars(res), vars(x)[2:4])
    expect_equal(specs(res), specs(x)[2:4, 2:4])

    res <- sel(x, 5, 4)
    expect_equal(as.matrix(res), orig[5, 4, drop=FALSE])
    expect_equal(vars(res), vars(x)[4])
    expect_equal(specs(res), specs(x)[4, 4, drop=FALSE])
})
