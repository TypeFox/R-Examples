context('sc_zip')

## dummy init list
dil <- list('dfvars' = TRUE,
            'select' = NULL,
            'filter' = NULL,
            'zip' = NULL,
            'year' = 2013)

test_that('Bad zip range, type, or missing', {
    expect_error(sc_zip(dil, 0000),
                 'Must provide a 5-digit zip code.')
    expect_error(sc_zip(dil, 000000),
                 'Must provide a 5-digit zip code.')
    expect_error(sc_zip(dil, '00000'),
                 'Must provide a 5-digit zip code.')
    expect_error(sc_zip(dil),
                 'Must provide a 5-digit zip code.')
})
