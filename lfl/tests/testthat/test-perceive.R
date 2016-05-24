test_that('perceive', {
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

    runPerceive <- function(...) {
        return(perceive(list(...), .vars, .specs))
    }

    expect_equal(runPerceive(c('Sm.d', 'Sm.a', 'Bi.c'),
                             c('Sm.d', 'VeSm.a', 'Bi.c'),
                             c('Sm.d', 'Sm.b', 'Sm.c')),
                 list(c('Sm.d', 'VeSm.a', 'Bi.c'),
                      c('Sm.d', 'Sm.b', 'Sm.c')))

    expect_equal(runPerceive(c('Sm.d', 'Sm.a', 'Bi.c'),
                             c('Sm.d', 'VeSm.a', 'Sm.c'),
                             c('Sm.d', 'Sm.b', 'Sm.c')),
                 list(c('Sm.d', 'Sm.a', 'Bi.c'),
                      c('Sm.d', 'VeSm.a', 'Sm.c'),
                      c('Sm.d', 'Sm.b', 'Sm.c')))

    expect_equal(runPerceive(c('Sm.d', 'Sm.a', 'Bi.c'),
                             c('Bi.d', 'VeSm.a', 'Bi.c'),
                             c('Sm.d', 'Sm.b', 'Sm.c')),
                 list(c('Bi.d', 'VeSm.a', 'Bi.c'),
                      c('Sm.d', 'Sm.b', 'Sm.c')))
})
