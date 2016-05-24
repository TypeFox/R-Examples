
###
### BASIC FUNCTIONALITY
###

context('Basic functionality of AppendCheckDigit')

test_that('AppendCheckDigit calls the appropriate function or throws errors', {
    expect_that(AppendCheckDigit(), throws_error())
    expect_that(AppendCheckDigit('1234'), throws_error())
    expect_that(AppendCheckDigit(method='Verhoeff'), throws_error())
    expect_that(AppendCheckDigit(1234, 'Verhoeff'), throws_error())
    expect_that(AppendCheckDigit('1234', 'foo'), throws_error())
    expect_that(AppendCheckDigit('1234', 'Verhoeff'), is_equivalent_to('12340'))
})



###
### VERHOEFF ALGORITHM
###

context('Verhoeff algorithm')

test_that('AppendCheckDigit.Verhoeff gives a warning when stripping non-digit characters', {
    expect_that(AppendCheckDigit('867-5309', 'Verhoeff'), gives_warning('Non-digit characters are disregarded in check digit calculation'))
})

test_that('AppendCheckDigit.Verhoeff returns correct values', {

    ## Missing values should receive no check digit
    expect_that(AppendCheckDigit('', 'Verhoeff'), is_equivalent_to(''))
    expect_that(AppendCheckDigit(as.character(NA), 'Verhoeff'), is_equivalent_to(as.character(NA)))
    expect_that(suppressWarnings(AppendCheckDigit(as.character(NaN), 'Verhoeff')), is_equivalent_to(as.character(NaN)))

    ## Character arguments
    expect_that(AppendCheckDigit('1', 'Verhoeff'), is_equivalent_to('15'))
    expect_that(AppendCheckDigit('1234', 'Verhoeff'), is_equivalent_to('12340'))
    expect_that(AppendCheckDigit('8675309', 'Verhoeff'), is_equivalent_to('86753098'))
    expect_that(AppendCheckDigit('9223372036854775808', 'Verhoeff'), is_equivalent_to('92233720368547758088'))

    ## Vectorized arguments
    expect_that(
        AppendCheckDigit(
            c('', '1', '1234', '8675309', '9223372036854775808'),
            'Verhoeff'
        ),
        is_equivalent_to(
            c('', '15', '12340', '86753098', '92233720368547758088')
        )
    )
})
