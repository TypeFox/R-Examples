context('sc_filter')

## dummy init list
dil <- list('dfvars' = TRUE,
            'select' = NULL,
            'filter' = NULL,
            'zip' = NULL,
            'year' = 2013)

test_that('Errors for non-init()', {
    expect_error(sc_filter(unitid == 99999),
                 'Chain not properly initialized. Be sure to start with sc_init().')
})

test_that('Errors for bad symbols', {
    expect_error(sc_filter(dil, unitid > 99999),
                 'Must use either \"==\" or \"!=\" in sc_filter.')
    expect_error(sc_filter(dil, unitid < 99999),
                 'Must use either \"==\" or \"!=\" in sc_filter.')
    expect_error(sc_filter(dil, unitid = 99999),
                 'Must use either \"==\" or \"!=\" in sc_filter.')
    expect_error(sc_filter(dil, unitid >= 99999),
                 'Must use either \"==\" or \"!=\" in sc_filter.')
    expect_error(sc_filter(dil, unitid <= 99999),
                 'Must use either \"==\" or \"!=\" in sc_filter.')
})


test_that('Error for bad variable names', {
    expect_error(sc_filter(dil, uniti == 99999))
})
