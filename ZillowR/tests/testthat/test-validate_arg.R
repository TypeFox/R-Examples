
context('validate_arg')

test_that("'required' validators accumulate errors as expected", {
    present_arg <- 'Here I am!'
    absent_arg <- NULL

    expect_identical(validate_arg(present_arg, required = FALSE), c())
    expect_identical(validate_arg(present_arg, required = TRUE), c())
    expect_identical(validate_arg(absent_arg, required = FALSE), c())
    expect_identical(validate_arg(absent_arg, required = TRUE), c("'absent_arg' is required"))
})

test_that("'class' validators accumulate errors as expected", {
    absent_arg <- NULL
    character_arg <- letters
    numeric_arg <- pi

    expect_identical(validate_arg(absent_arg, class = 'character'), c())
    expect_identical(validate_arg(absent_arg, class = 'numeric'), c())
    expect_identical(validate_arg(character_arg, class = 'character'), c())
    expect_identical(validate_arg(character_arg, class = 'numeric'), c("'character_arg' must be numeric"))
    expect_identical(validate_arg(numeric_arg, class = 'character'), c("'numeric_arg' must be character"))
    expect_identical(validate_arg(numeric_arg, class = 'numeric'), c())
})

test_that("'length_min' validators accumulate errors as expected", {
    absent_arg <- NULL
    len1_arg <- 1
    len2_arg <- 1:2

    expect_identical(validate_arg(absent_arg, length_min = 1), c())
    expect_identical(validate_arg(absent_arg, length_min = 2), c())
    expect_identical(validate_arg(len1_arg, length_min = 1), c())
    expect_identical(validate_arg(len1_arg, length_min = 2), c("'len1_arg' must be at least length 2"))
    expect_identical(validate_arg(len2_arg, length_min = 1), c())
    expect_identical(validate_arg(len2_arg, length_min = 2), c())
})

test_that("'length_max' validators accumulate errors as expected", {
    absent_arg <- NULL
    len1_arg <- 1
    len2_arg <- 1:2

    expect_identical(validate_arg(absent_arg, length_max = 1), c())
    expect_identical(validate_arg(absent_arg, length_max = 2), c())
    expect_identical(validate_arg(len1_arg, length_max = 1), c())
    expect_identical(validate_arg(len1_arg, length_max = 2), c())
    expect_identical(validate_arg(len2_arg, length_max = 1), c("'len2_arg' must be no more than length 1"))
    expect_identical(validate_arg(len2_arg, length_max = 2), c())
})

test_that("'inclusion' validators accumulate errors as expected", {
    absent_arg <- NULL
    letters_01 <- head(letters, 1)
    letters_26 <- letters
    numbers_01 <- 1
    numbers_26 <- 1:26
    letters_and_numbers <- c(letters_26, numbers_26)

    expect_identical(validate_arg(absent_arg, inclusion = letters), c())
    expect_identical(validate_arg(absent_arg, inclusion = 1:26), c())
    expect_identical(validate_arg(letters_01, inclusion = letters), c())
    expect_identical(validate_arg(letters_01, inclusion = 1:26), c(sprintf("'letters_01' must be one of: %s", paste(1:26, collapse = ', '))))
    expect_identical(validate_arg(letters_26, inclusion = letters), c())
    expect_identical(validate_arg(letters_26, inclusion = 1:26), c(sprintf("'letters_26' must be one of: %s", paste(1:26, collapse = ', '))))
    expect_identical(validate_arg(numbers_01, inclusion = letters), c(sprintf("'numbers_01' must be one of: %s", paste(letters, collapse = ', '))))
    expect_identical(validate_arg(numbers_01, inclusion = 1:26), c())
    expect_identical(validate_arg(numbers_26, inclusion = letters), c(sprintf("'numbers_26' must be one of: %s", paste(letters, collapse = ', '))))
    expect_identical(validate_arg(numbers_26, inclusion = 1:26), c())
    expect_identical(validate_arg(letters_and_numbers, inclusion = letters), c(sprintf("'letters_and_numbers' must be one of: %s", paste(letters, collapse = ', '))))
    expect_identical(validate_arg(letters_and_numbers, inclusion = 1:26), c(sprintf("'letters_and_numbers' must be one of: %s", paste(1:26, collapse = ', '))))
})

test_that("'exclusion' validators accumulate errors as expected", {
    absent_arg <- NULL
    letters_01 <- head(letters, 1)
    letters_26 <- letters
    numbers_01 <- 1
    numbers_26 <- 1:26
    letters_and_numbers <- c(letters_26, numbers_26)

    expect_identical(validate_arg(absent_arg, exclusion = letters), c())
    expect_identical(validate_arg(absent_arg, exclusion = 1:26), c())
    expect_identical(validate_arg(letters_01, exclusion = letters), c(sprintf("'letters_01' must not be any of: %s", paste(letters, collapse = ', '))))
    expect_identical(validate_arg(letters_01, exclusion = 1:26), c())
    expect_identical(validate_arg(letters_26, exclusion = letters), c(sprintf("'letters_26' must not be any of: %s", paste(letters, collapse = ', '))))
    expect_identical(validate_arg(letters_26, exclusion = 1:26), c())
    expect_identical(validate_arg(numbers_01, exclusion = letters), c())
    expect_identical(validate_arg(numbers_01, exclusion = 1:26), c(sprintf("'numbers_01' must not be any of: %s", paste(1:26, collapse = ', '))))
    expect_identical(validate_arg(numbers_26, exclusion = letters), c())
    expect_identical(validate_arg(numbers_26, exclusion = 1:26), c(sprintf("'numbers_26' must not be any of: %s", paste(1:26, collapse = ', '))))
    expect_identical(validate_arg(letters_and_numbers, exclusion = letters), c(sprintf("'letters_and_numbers' must not be any of: %s", paste(letters, collapse = ', '))))
    expect_identical(validate_arg(letters_and_numbers, exclusion = 1:26), c(sprintf("'letters_and_numbers' must not be any of: %s", paste(1:26, collapse = ', '))))
})

test_that("'format' validators accumulate errors as expected", {
    absent_arg <- NULL
    beanstalk_arg <- c('fee', 'fi', 'fo', 'fum')
    numeric_arg <- 1:3

    expect_identical(validate_arg(absent_arg, format = '^f'), c())
    expect_identical(validate_arg(absent_arg, format = '^f[aeiou]$'), c())
    expect_identical(validate_arg(absent_arg, format = '^\\d+$'), c())
    expect_identical(validate_arg(beanstalk_arg, format = '^f'), c())
    expect_identical(validate_arg(beanstalk_arg, format = '^f[aeiou]$'), c("'beanstalk_arg' must match: ^f[aeiou]$"))
    expect_identical(validate_arg(beanstalk_arg, format = '^\\d+$'), c("'beanstalk_arg' must match: ^\\d+$"))
    expect_identical(validate_arg(numeric_arg, format = '^f'), c("'numeric_arg' must match: ^f"))
    expect_identical(validate_arg(numeric_arg, format = '^f[aeiou]$'), c("'numeric_arg' must match: ^f[aeiou]$"))
    expect_identical(validate_arg(numeric_arg, format = '^\\d+$'), c())
})

test_that("'value_min' validators accumulate errors as expected", {
    absent_arg <- NULL
    value1_arg <- 1
    value2_arg <- 1:2

    expect_identical(validate_arg(absent_arg, value_min = 1), c())
    expect_identical(validate_arg(absent_arg, value_min = 2), c())
    expect_identical(validate_arg(value1_arg, value_min = 1), c())
    expect_identical(validate_arg(value1_arg, value_min = 2), c("'value1_arg' values must be greater than or equal to 2"))
    expect_identical(validate_arg(value2_arg, value_min = 1), c())
    expect_identical(validate_arg(value2_arg, value_min = 2), c("'value2_arg' values must be greater than or equal to 2"))
})

test_that("'value_max' validators accumulate errors as expected", {
    absent_arg <- NULL
    value1_arg <- 1
    value2_arg <- 1:2

    expect_identical(validate_arg(absent_arg, value_max = 1), c())
    expect_identical(validate_arg(absent_arg, value_max = 2), c())
    expect_identical(validate_arg(value1_arg, value_max = 1), c())
    expect_identical(validate_arg(value1_arg, value_max = 2), c())
    expect_identical(validate_arg(value2_arg, value_max = 1), c("'value2_arg' values must be less than or equal to 1"))
    expect_identical(validate_arg(value2_arg, value_max = 2), c())
})
