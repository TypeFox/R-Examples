context('sc_get')

test_that('Errors for non-init()', {
    expect_error(sc_get(),
                 'Chain not properly initialized. Be sure to start with sc_init().')
})

