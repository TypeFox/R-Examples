test_that(
  "test.character_to_list_of_integer_vectors.strings.returns_list_of_integer_vectors", 
  {
    x <- c("12345", "1b3d5", "abcde", NA, "", " ", " 2 4 ")
    expected <- list(1:5, c(1L, NA, 3L, NA, 5L), rep.int(NA_integer_, 5L), 
                     NA_integer_, integer(), NA_integer_, c(NA, 2L, NA, 4L, NA))
    names(expected) <- x
    expect_identical(suppressWarnings(character_to_list_of_integer_vectors(x)), 
                     expected)
    expect_warning(character_to_list_of_integer_vectors(x))
  })

test_that(
  "test.d.single_input.returns_valid_regex",
  {
    lo <- 1:5
    expected <- paste0(
      "[[:digit:]]", 
      c("", paste0("{", lo[-1], "}"))
    )
    actual <- d(lo)
    expect_identical(actual, expected)
  }
)

test_that(
  "test.d.both_inputs.returns_valid_regex",
  {
    lo <- 1:5
    hi <- 6:8
    expected <- paste0("[[:digit:]]{", lo, ",", hi, "}")
    actual <- d(lo, hi)
    expect_identical(actual, expected)
  }
)

test_that(
  "test.d.infinite_hi.returns_valid_regex",
  {
    lo <- 0:2
    hi <- Inf
    expected <- paste0("[[:digit:]]", c("*", "+", "{2,}"))
    actual <- d(lo, hi)
    expect_identical(actual, expected)
  }
)

test_that(
  "test.matches_regex.strings.returns_true_when_string_matches_regex", 
  {
    rx <- "foo"
    x <- c("foo", "fooo", "fo", "", "FOO", NA)
    expected <- c(TRUE, TRUE, FALSE, FALSE, TRUE, NA)
    names(expected) <- x
    expect_equal(matches_regex(x, rx), expected)
  }) 
