test_that('fcut of numeric by single triangle', {
    x <- 0:100
    res <- fcut(x, 
                 breaks=c(0, 50, 100),
                 type='triangle')

    expect_true(is.matrix(res))
    expect_equal(ncol(res), 1)
    expect_equal(nrow(res), 101)
    expect_equal(colnames(res), 'x.1')
    expect_true(inherits(res, 'fsets'))
    expect_equivalent(vars(res), 'x')
    expect_equal(names(vars(res)), colnames(res))
    expect_equal(specs(res), matrix(0, 
                                    nrow=1,
                                    ncol=1,
                                    dimnames=list(colnames(res), colnames(res))))

    expect_equivalent(res[1, 1], 0)
    expect_equivalent(res[26, 1], 0.5)
    expect_equivalent(res[51, 1], 1)
    expect_equivalent(res[76, 1], 0.5)
    expect_equivalent(res[101, 1], 0)
    expect_true(is.fsets(res))
})


test_that('fcut of numeric by multiple triangles', {
    x <- 0:100
    res <- fcut(x, 
                 breaks=c(0, 25, 50, 75, 100),
                 type='triangle')

    expect_true(is.matrix(res))
    expect_equal(ncol(res), 3)
    expect_equal(nrow(res), 101)
    expect_equal(colnames(res), c('x.1', 'x.2', 'x.3'))
    expect_true(inherits(res, 'fsets'))
    expect_equivalent(vars(res), rep('x', 3))
    expect_equal(names(vars(res)), colnames(res))
    expect_equal(specs(res), matrix(rep(0, 3*3), 
                                    nrow=3,
                                    ncol=3,
                                    dimnames=list(colnames(res), colnames(res))))

    expect_equivalent(res[1, 1], 0)
    expect_equivalent(res[26, 1], 1)
    expect_equivalent(res[51, 1], 0)
    expect_equivalent(res[76, 1], 0)
    expect_equivalent(res[101, 1], 0)

    expect_equivalent(res[1, 2], 0)
    expect_equivalent(res[26, 2], 0)
    expect_equivalent(res[51, 2], 1)
    expect_equivalent(res[76, 2], 0)
    expect_equivalent(res[101, 2], 0)

    expect_equivalent(res[1, 3], 0)
    expect_equivalent(res[26, 3], 0)
    expect_equivalent(res[51, 3], 0)
    expect_equivalent(res[76, 3], 1)
    expect_equivalent(res[101, 3], 0)
    expect_true(is.fsets(res))
})


test_that('fcut of numeric with merge 1:3', {
    x <- 0:100
    res <- fcut(x, 
                breaks=c(0, 25, 50, 75, 100),
                type='triangle',
                merge=1:3)
   
    expect_true(is.matrix(res))
    expect_equal(ncol(res), 6)
    expect_equal(nrow(res), 101)
    expect_equal(colnames(res), c('x.1', 'x.2', 'x.3',
                                  'x.1|x.2', 'x.2|x.3',
                                  'x.1|x.2|x.3'))
    expect_true(inherits(res, 'fsets'))
    expect_equivalent(vars(res), rep('x', 6))
    expect_equal(names(vars(res)), colnames(res))
    expect_equal(specs(res), matrix(c(0,0,0,1,0,1,
                                      0,0,0,1,1,1,
                                      0,0,0,0,1,1,
                                      0,0,0,0,0,1,
                                      0,0,0,0,0,1,
                                      0,0,0,0,0,0),
                                    byrow=TRUE,
                                    nrow=6,
                                    dimnames=list(colnames(res), colnames(res))))
    expect_true(is.fsets(res))
})


test_that('fcut of numeric with merge 2', {
    x <- 0:100
    res <- fcut(x, 
                breaks=c(0, 25, 50, 75, 100),
                type='triangle',
                merge=2)
   
    expect_true(is.matrix(res))
    expect_equal(ncol(res), 2)
    expect_equal(nrow(res), 101)
    expect_equal(colnames(res), c('x.1|x.2', 'x.2|x.3'))
    expect_true(inherits(res, 'fsets'))
    expect_equivalent(vars(res), rep('x', 2))
    expect_equal(names(vars(res)), colnames(res))
    expect_equal(specs(res), matrix(rep(0, 2*2), 
                                    nrow=2,
                                    ncol=2,
                                    dimnames=list(colnames(res), colnames(res))))

    expect_equivalent(res[1, 1], 0)
    expect_equivalent(res[26, 1], 1)
    expect_equivalent(res[30, 1], 1)
    expect_equivalent(res[51, 1], 1)
    expect_equivalent(res[76, 1], 0)
    expect_equivalent(res[101, 1], 0)

    expect_equivalent(res[1, 2], 0)
    expect_equivalent(res[26, 2], 0)
    expect_equivalent(res[51, 2], 1)
    expect_equivalent(res[60, 2], 1)
    expect_equivalent(res[76, 2], 1)
    expect_equivalent(res[101, 2], 0)
    expect_true(is.fsets(res))
})


test_that('fcut of numeric with merge 1,3', {
    x <- 0:100
    res <- fcut(x, 
                breaks=c(0, 25, 50, 75, 100),
                type='triangle',
                merge=c(1,3))
   
    expect_true(is.matrix(res))
    expect_equal(ncol(res), 4)
    expect_equal(nrow(res), 101)
    expect_equal(colnames(res), c('x.1', 'x.2', 'x.3',
                                  'x.1|x.2|x.3'))
    expect_true(inherits(res, 'fsets'))
    expect_equivalent(vars(res), rep('x', 4))
    expect_equal(names(vars(res)), colnames(res))
    expect_equal(specs(res), matrix(c(0,0,0,1,
                                      0,0,0,1,
                                      0,0,0,1,
                                      0,0,0,0),
                                    byrow=TRUE,
                                    nrow=4,
                                    dimnames=list(colnames(res), colnames(res))))
    expect_true(is.fsets(res))
})
