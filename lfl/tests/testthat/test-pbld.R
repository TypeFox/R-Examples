test_that('pbld', {
    # init fsets
    .vars <- c(rep('b', 3),
               rep('c', 3))
    names(.vars) <- paste(rep(c('VeSm', 'Sm', 'Bi'), times=2),
                          rep(c('b', 'c'), each=3),
                          sep='.')


    .specs <- matrix(c(0,1,0, 0,0,0,
                       0,0,0, 0,0,0,
                       0,0,0, 0,0,0,
 
                       0,0,0, 0,1,0,
                       0,0,0, 0,0,0,
                       0,0,0, 0,0,0),
                     byrow=TRUE,
                     ncol=6)
    colnames(.specs) <- names(.vars)
    rownames(.specs) <- names(.vars)

    x <- matrix(runif(18), ncol=6, nrow=3)
    colnames(x) <- names(.vars)

    x <- fsets(x, vars=.vars, specs=.specs)

    # init rules
    rules <- list(c('Sm.b', 'VeSm.c'),
                  c('Sm.b', 'Sm.c'),
                  c('Bi.b', 'Bi.b'),
                  c('Bi.b', 'Sm.c', 'Sm.b'))

    # init values
    values <- 0:10 / 10

    # init partition
    partition <- lcut3(data.frame(b=values))


    res <- pbld(x, rules, partition, values)
    expect_equal(nrow(x), length(res))
    #expect_true(FALSE)
})


test_that('pbld with empty rulebase', {
    # init fsets
    .vars <- c(rep('b', 3),
               rep('c', 3))
    names(.vars) <- paste(rep(c('VeSm', 'Sm', 'Bi'), times=2),
                          rep(c('b', 'c'), each=3),
                          sep='.')


    .specs <- matrix(c(0,1,0, 0,0,0,
                       0,0,0, 0,0,0,
                       0,0,0, 0,0,0,
 
                       0,0,0, 0,1,0,
                       0,0,0, 0,0,0,
                       0,0,0, 0,0,0),
                     byrow=TRUE,
                     ncol=6)
    colnames(.specs) <- names(.vars)
    rownames(.specs) <- names(.vars)

    x <- matrix(runif(18), ncol=6, nrow=3)
    colnames(x) <- names(.vars)

    x <- fsets(x, vars=.vars, specs=.specs)

    # init rules
    rules <- list()

    # init values
    values <- 0:10 / 10

    # init partition
    partition <- lcut3(data.frame(b=values))


    res <- pbld(x, rules, partition, values)
    expect_equal(nrow(x), length(res))
})
