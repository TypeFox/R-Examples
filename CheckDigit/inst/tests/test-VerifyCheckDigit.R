
###
### SETUP
###

SingleDigitSubstitutions <- function(x) {
    stopifnot(is.character(x) & length(x) == 1 & all(grepl('^\\d*$', x)))

    y <- c()

    digits <- as.integer(unlist(strsplit(x, '')))
    for (i in seq(digits)) {
        for (j in setdiff(0:9, digits[i])) {
            new_digits <- digits
            new_digits[i] <- j
            y <- c(y, paste(new_digits, collapse=''))
            rm(new_digits)
        }
        rm(j)
    }
    rm(i)

    return(y)
}

SingleAdjacentTranspositions <- function(x) {
    stopifnot(is.character(x) & length(x) == 1 & all(grepl('^\\d*$', x)))

    y <- c()

    if (nchar(x) > 1) {
        digits <- as.integer(unlist(strsplit(x, '')))
        for (i in seq(2, length(digits))) {
            index <- seq(digits)
            index[i-1] <- index[i-1] + 1
            index[i] <- index[i] - 1
            y <- c(y, paste(digits[index], collapse=''))
            rm(index)
        }
        rm(i)
    }

    return(y)
}



###
### BASIC FUNCTIONALITY
###

context('Basic functionality of VerifyCheckDigit')

test_that('VerifyCheckDigit calls the appropriate function or throws errors', {
    expect_that(VerifyCheckDigit(), throws_error())
    expect_that(VerifyCheckDigit('12340'), throws_error())
    expect_that(VerifyCheckDigit(method='Verhoeff'), throws_error())
    expect_that(VerifyCheckDigit(12340, 'Verhoeff'), throws_error())
    expect_that(VerifyCheckDigit('12340', 'foo'), throws_error())
    expect_that(VerifyCheckDigit('12340', 'Verhoeff'), is_equivalent_to(TRUE))
})



###
### VERHOEFF ALGORITHM
###

context('Verhoeff algorithm')

test_that('VerifyCheckDigit.Verhoeff gives a warning when stripping non-digit characters', {
    expect_that(VerifyCheckDigit('867-5309', 'Verhoeff'), gives_warning('Non-digit characters are disregarded in check digit calculation'))
})

test_that('VerifyCheckDigit.Verhoeff returns correct responses', {

    ## Missing values should return FALSE
    expect_that(VerifyCheckDigit('', 'Verhoeff'), is_equivalent_to(FALSE))
    expect_that(VerifyCheckDigit(as.character(NA), 'Verhoeff'), is_equivalent_to(FALSE))
    expect_that(suppressWarnings(VerifyCheckDigit(as.character(NaN), 'Verhoeff')), is_equivalent_to(FALSE))

    ## Character arguments
    expect_that(VerifyCheckDigit('15', 'Verhoeff'), is_equivalent_to(TRUE))
    expect_that(VerifyCheckDigit('12340', 'Verhoeff'), is_equivalent_to(TRUE))
    expect_that(VerifyCheckDigit('86753098', 'Verhoeff'), is_equivalent_to(TRUE))
    expect_that(VerifyCheckDigit('92233720368547758088', 'Verhoeff'), is_equivalent_to(TRUE))

    ## Vectorized arguments
    expect_that(
        VerifyCheckDigit(
            c('', '15', '12340', '86753098', '92233720368547758088'),
            'Verhoeff'
        ),
        is_equivalent_to(
            c(FALSE, TRUE, TRUE, TRUE, TRUE)
        )
    )
})

test_that('Verhoeff check digit detects all single digit substitutions', {
    for (i in SingleDigitSubstitutions('86753098')) {
        expect_that(VerifyCheckDigit(i, 'Verhoeff'), is_equivalent_to(FALSE))
    }
    rm(i)
})

test_that('Verhoeff check digit detects all single adjacent transpositions', {
    for (i in SingleAdjacentTranspositions('86753098')) {
        expect_that(VerifyCheckDigit(i, 'Verhoeff'), is_equivalent_to(FALSE))
    }
    rm(i)
})
