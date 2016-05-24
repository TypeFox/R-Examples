context("h")

test_that_for_each_copula("h() in [0,1]", {
    u <- h(copula, XV[ , 1], XV[ , 2])
    expect_true(all(is.finite(u)))
    expect_true(all(u >= 0 & u <= 1))
})

test_that_for_each_copula("h() is correct", {
    skip_on_cran()

    XV <- subset(XV, X >= 0.25 & X <= 0.75 & V >= 0.25 & V <= 0.75)
    u <- h(copula, XV[ , 1], XV[ , 2])
    uu <- vines:::hCopula(copula, XV[ , 1], XV[ , 2],
                          eps = .Machine$double.eps^0.5)
    expect_equal(u, uu, tolerance = tol)
})

test_that_for_each_copula("h(0, v) == 0", {
    u <- h(copula, rep(0, n), V)
    expect_equal(u, rep(0, n), tolerance = tol)
})

test_that_for_each_copula("h(1, v) == 1", {
    u <- h(copula, rep(1, n), V)
    expect_equal(u, rep(1, n), tolerance = tol)
})
