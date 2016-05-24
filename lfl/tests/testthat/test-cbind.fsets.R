test_that('cbind of two fsets', {
    x <- 0:10 * 10
    y <- 0:10 * 10
    a <- fcut(x, breaks=c(0, 25, 50, 75, 100), merge=1:2)    #x1, x2, x3, x1|2, x2|3
    b <- fcut(y, breaks=c(0, 30, 70, 80, 100), merge=c(1,3)) #y1, y2, y3, y1|2|3
    res <- cbind.fsets(a, b)


    expectedNames <- c('x.1', 'x.2', 'x.3',
                       'x.1|x.2', 'x.2|x.3',
                       'y.1', 'y.2', 'y.3',
                       'y.1|y.2|y.3')

    expect_true(is.matrix(res))
    expect_equal(ncol(res), 9)
    expect_equal(nrow(res), 11)
    expect_equal(colnames(res), expectedNames)
    expect_true(inherits(res, 'fsets'))

    expectedVars <- c(rep('x', 5), rep('y', 4))
    names(expectedVars) <- expectedNames

    expect_equal(vars(res), expectedVars)
    expect_equal(specs(res), matrix(c(0,0,0,1,0, 0,0,0,0,
                                      0,0,0,1,1, 0,0,0,0,
                                      0,0,0,0,1, 0,0,0,0,
                                      0,0,0,0,0, 0,0,0,0,
                                      0,0,0,0,0, 0,0,0,0,

                                      0,0,0,0,0, 0,0,0,1,
                                      0,0,0,0,0, 0,0,0,1,
                                      0,0,0,0,0, 0,0,0,1,
                                      0,0,0,0,0, 0,0,0,0),
                                    byrow=TRUE,
                                    nrow=9,
                                    dimnames=list(expectedNames, expectedNames)))
    expect_equal(res[, 1:5], a[, 1:ncol(a)])  # 'a[, ...]' removes all attributes from 'a'
    expect_equal(res[, 6:9], b[, 1:ncol(b)])  # 'b[, ...]' removes all attributes from 'b'
    expect_true(is.fsets(res))
})



test_that('cbind of two fsets and NULL', {
    x <- 0:10 * 10
    y <- 0:10 * 10
    a <- fcut(x, breaks=c(0, 25, 50, 75, 100), merge=1:2)    #x1, x2, x3, x1|2, x2|3
    b <- fcut(y, breaks=c(0, 30, 70, 80, 100), merge=c(1,3)) #y1, y2, y3, y1|2|3
    res <- cbind.fsets(a, b, NULL)

    expectedNames <- c('x.1', 'x.2', 'x.3',
                       'x.1|x.2', 'x.2|x.3',
                       'y.1', 'y.2', 'y.3',
                       'y.1|y.2|y.3')

    expect_true(is.matrix(res))
    expect_equal(ncol(res), 9)
    expect_equal(nrow(res), 11)
    expect_equal(colnames(res), expectedNames)
    expect_true(inherits(res, 'fsets'))

    expectedVars <- c(rep('x', 5), rep('y', 4))
    names(expectedVars) <- expectedNames

    expect_equal(vars(res), expectedVars)
    expect_equal(specs(res), matrix(c(0,0,0,1,0, 0,0,0,0,
                                      0,0,0,1,1, 0,0,0,0,
                                      0,0,0,0,1, 0,0,0,0,
                                      0,0,0,0,0, 0,0,0,0,
                                      0,0,0,0,0, 0,0,0,0,

                                      0,0,0,0,0, 0,0,0,1,
                                      0,0,0,0,0, 0,0,0,1,
                                      0,0,0,0,0, 0,0,0,1,
                                      0,0,0,0,0, 0,0,0,0),
                                    byrow=TRUE,
                                    nrow=9,
                                    dimnames=list(expectedNames, expectedNames)))
    expect_equal(res[, 1:5], a[, 1:ncol(a)])  # 'a[, ...]' removes all attributes from 'a'
    expect_equal(res[, 6:9], b[, 1:ncol(b)])  # 'b[, ...]' removes all attributes from 'b'
    expect_true(is.fsets(res))
})

