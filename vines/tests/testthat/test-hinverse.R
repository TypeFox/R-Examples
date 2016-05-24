context("hinverse")

test_that_for_each_copula("hinverse() in [0,1]", {
    x <- hinverse(copula, XV[ , 1], XV[ , 2])
    expect_true(all(is.finite(x)))
    expect_true(all(x >= 0 & x <= 1))
})

test_that_for_each_copula("hinverse() is correct", {
    skip_on_cran()

    XV <- subset(XV, X >= 0.25 & X <= 0.75 & V >= 0.25 & V <= 0.75)
    x <- hinverse(copula, XV[ , 1], XV[ , 2])
    xx <- vines:::hinverseCopula(copula, XV[ , 1], XV[ , 2],
                                 eps = .Machine$double.eps^0.5)
    expect_equal(x, xx, tolerance = tol)
})

test_that_for_each_copula("hinverse(0, v) == 0", {
    x <- hinverse(copula, rep(0, n), V)
    expect_equal(x, rep(0, n), tolerance = tol)
})

test_that_for_each_copula("hinverse(1, v) == 1", {
    x <- hinverse(copula, rep(1, n), V)
    expect_equal(x, rep(1, n), tolerance = tol)
})
