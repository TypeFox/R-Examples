test_that('is.specific', {
    .vars <- c(rep('a', 3),
              rep('b', 3),
              rep('c', 3),
              rep('d', 3))
    names(.vars) <- paste(rep(c('VeSm', 'Sm', 'Bi'), times=4),
                         rep(c('a', 'b', 'c', 'd'), each=3),
                         sep='.')

    .specs <- matrix(c(0,1,0, 0,0,0, 0,0,0, 0,0,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0,

                      0,0,0, 0,1,0, 0,0,0, 0,0,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0,

                      0,0,0, 0,0,0, 0,1,0, 0,0,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0,

                      0,0,0, 0,0,0, 0,0,0, 0,1,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0,
                      0,0,0, 0,0,0, 0,0,0, 0,0,0),
                    byrow=TRUE,
                    ncol=12)
    colnames(.specs) = names(.vars)
    rownames(.specs) = names(.vars)

    # equal rules
    expect_true(is.specific(c('VeSm.a', 'Bi.c'),
                            c('VeSm.a', 'Bi.c'),
                            .vars, .specs))

    # the same but single item that is more specific
    expect_true(is.specific(c('VeSm.a', 'Bi.c', 'Sm.d'),
                            c('Sm.a', 'Bi.c', 'Sm.d'),
                            .vars, .specs))

    # "any" is less specific
    expect_true(is.specific(c('VeSm.a', 'Bi.c', 'Sm.d'),
                            c('VeSm.a', 'Bi.c'),
                            .vars, .specs))

    # "any" + other
    expect_true(is.specific(c('VeSm.a', 'Bi.c', 'Sm.d'),
                            c('Sm.a', 'Bi.c'),
                            .vars, .specs))

    # everything is more specific than empty rule
    expect_true(is.specific(c('VeSm.a', 'Bi.c', 'Sm.d'),
                            NULL,
                            .vars, .specs))

    # null rules are <=
    expect_true(is.specific(NULL,
                            NULL,
                            .vars, .specs))

    # the same but single item that is more specific
    expect_false(is.specific(c('Sm.a', 'Bi.c', 'Sm.d'),
                             c('VeSm.a', 'Bi.c', 'Sm.d'),
                             .vars, .specs))

    # "any" is less specific
    expect_false(is.specific(c('VeSm.a', 'Bi.c'),
                             c('VeSm.a', 'Bi.c', 'Sm.d'),
                             .vars, .specs))

    # "any" + other
    expect_false(is.specific(c('Sm.a', 'Bi.c'),
                             c('VeSm.a', 'Bi.c', 'Sm.d'),
                             .vars, .specs))

    # everything is more specific than empty rule
    expect_false(is.specific(NULL,
                             c('VeSm.a', 'Bi.c', 'Sm.d'),
                             .vars, .specs))

    # different .vars are incomparable
    expect_false(is.specific(c('Sm.a'),
                             c('Bi.c'),
                             .vars, .specs))
    expect_false(is.specific(c('Bi.c'),
                             c('Sm.a'),
                             .vars, .specs))



    expect_false(is.specific(c('VeSm.a', 'Sm.c'),
                             c('Sm.a', 'Bi.c'),
                             .vars, .specs))
    expect_false(is.specific(c('Sm.b', 'Sm.d'),
                             c('Sm.a', 'Bi.c'),
                             .vars, .specs))
})
