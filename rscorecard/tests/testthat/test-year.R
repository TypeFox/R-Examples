context('sc_year')

## dummy init list
dil <- list('dfvars' = TRUE,
            'select' = NULL,
            'filter' = NULL,
            'zip' = NULL,
            'year' = 2013)

test_that('Errors for non-init()', {
    expect_error(sc_year(2000),
                 'Chain not properly initialized. Be sure to start with sc_init().')
})

test_that('Bad year range, type, or missing', {
    expect_error(sc_year(dil, 1847),
                 'Must provide a 4-digit year in 1900s or 2000s.')
    expect_error(sc_year(dil, 2100),
                 'Must provide a 4-digit year in 1900s or 2000s.')
    expect_error(sc_year(dil, '2000'),
                 'Must provide a 4-digit year in 1900s or 2000s.')
    expect_error(sc_year(dil),
                 'Must provide a 4-digit year in 1900s or 2000s.')

})
