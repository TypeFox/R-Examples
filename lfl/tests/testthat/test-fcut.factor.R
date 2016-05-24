test_that('fcut of factor', {
    x <- as.factor(c('a', 'b', 'a', 'c', 'c', 'b', 'c'))

    res <- fcut(x)


    expect_true(is.matrix(res))
    expect_equal(ncol(res), 3)
    expect_equal(nrow(res), 7)
    expect_equal(colnames(res), c('x.a', 'x.b', 'x.c'))
    expect_true(inherits(res, 'fsets'))
    expect_equivalent(vars(res), rep('x', 3))
    expect_equal(names(vars(res)), colnames(res))
    expect_equal(specs(res), matrix(rep(0, 3*3), 
                                     nrow=3,
                                     ncol=3,
                                     dimnames=list(colnames(res), colnames(res))))
    expect_true(is.fsets(res))

    expect_equivalent(res[1, 1], 1)
    expect_equivalent(res[1, 2], 0)
    expect_equivalent(res[1, 3], 0)

    expect_equivalent(res[2, 1], 0)
    expect_equivalent(res[2, 2], 1)
    expect_equivalent(res[2, 3], 0)

    expect_equivalent(res[3, 1], 1)
    expect_equivalent(res[3, 2], 0)
    expect_equivalent(res[3, 3], 0)

    expect_equivalent(res[4, 1], 0)
    expect_equivalent(res[4, 2], 0)
    expect_equivalent(res[4, 3], 1)

    expect_equivalent(res[5, 1], 0)
    expect_equivalent(res[5, 2], 0)
    expect_equivalent(res[5, 3], 1)

    expect_equivalent(res[6, 1], 0)
    expect_equivalent(res[6, 2], 1)
    expect_equivalent(res[6, 3], 0)

    expect_equivalent(res[7, 1], 0)
    expect_equivalent(res[7, 2], 0)
    expect_equivalent(res[7, 3], 1)
})
