context("repeated, bad lo")

test_that(
  "repeated with lo < 0 throws an error",
  {
    expect_error(
      repeated("foo", -1),
      "lo has negative values"
    )
  }
)

test_that(
  "repeated with missing/NaN/infinite lo throws an error",
  {
    for(lo in c(NA, NaN, Inf))
    {
      expect_error(
        repeated("foo", lo),
        "lo has missing or infinite values",
        info = paste("lo =", lo)
      )
    }
  }
)

context("repeated, once")

test_that(
  "repeated with lo, hi missing, char_class = FALSE returns the input",
  {
    expected <- as.regex("foo")
    actual <- repeated("foo", char_class = FALSE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 1, hi missing, char_class = FALSE returns the input",
  {
    expected <- as.regex("foo")
    actual <- repeated("foo", 1, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

context("repeated, {} output")

test_that(
  "repeated with lo = 2, hi missing, char_class = FALSE returns the input and {2}",
  {
    expected <- as.regex("foo{2}")
    actual <- repeated("foo", 2, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 2, hi = 5, char_class = FALSE returns the input and {2,5}",
  {
    expected <- as.regex("foo{2,5}")
    actual <- repeated("foo", 2, 5, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 2, hi = Inf, char_class = FALSE returns the input and {2,}",
  {
    expected <- as.regex("foo{2,}")
    actual <- repeated("foo", 2, Inf, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 2, hi = 2, char_class = FALSE returns the input and {2}",
  {
    expected <- as.regex("foo{2}")
    actual <- repeated("foo", 2, 2, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

context("repeated, impossible no. of reps")

test_that(
  "repeated with hi < lo, throws an error",
  {
    expect_error(
      repeated("foo", 2, 1),
      "hi has values that are less than the corresponding values in lo"
    )
  }
)

context("repeated, special no. of reps, char_class = FALSE")

test_that(
  "repeated with lo = 1, hi = Inf, char_class = FALSE returns the input and +",
  {
    expected <- as.regex("foo+")
    actual <- repeated("foo", 1, Inf, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 0, hi = Inf, char_class = FALSE returns the input and *",
  {
    expected <- as.regex("foo+")
    actual <- repeated("foo", 1, Inf, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 0, hi = 1, char_class = FALSE returns the input and ?",
  {
    expected <- as.regex("foo?")
    actual <- repeated("foo", 0, 1, char_class = FALSE)
    expect_equal(actual, expected)
  }
)

context("repeated, char_class = TRUE")

test_that(
  "repeated with lo = 2, hi = 5, char_class = TRUE repeats + wraps input in []",
  {
    expected <- as.regex("[foo]{2,5}")
    actual <- repeated("foo", 2, 5, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

# Next 3 tests mostly checking that you don't get double class wrapping [[]]

test_that(
  "repeated with lo = 1, hi = Inf, char_class = TRUE repeats + wraps input in []",
  {
    expected <- as.regex("[foo]+")
    actual <- repeated("foo", 1, Inf, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 0, hi = Inf, char_class = TRUE repeats + wraps input in []",
  {
    expected <- as.regex("[foo]*")
    actual <- repeated("foo", 0, Inf, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 0, hi = 1, char_class = TRUE repeats + wraps input in []",
  {
    expected <- as.regex("[foo]?")
    actual <- repeated("foo", 0, 1, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

context("repeated, lazy arg")

test_that(
  "repeated with lo = 2, hi = 5, lazy = TRUE repeats + ?",
  {
    expected <- as.regex("[foo]{2,5}?")
    actual <- repeated("foo", 2, 5, lazy = TRUE, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 1, hi = Inf, lazy = TRUE repeats + ?",
  {
    expected <- as.regex("[foo]+?")
    actual <- repeated("foo", 1, Inf, lazy = TRUE, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 0, hi = Inf, lazy = TRUE repeats + ?",
  {
    expected <- as.regex("[foo]*?")
    actual <- repeated("foo", 0, Inf, lazy = TRUE, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "repeated with lo = 0, hi = 1, lazy = TRUE repeats + ?",
  {
    expected <- as.regex("[foo]??") # TODO: Is this right, or should it only be 1 '?'
    actual <- repeated("foo", 0, 1, lazy = TRUE, char_class = TRUE)
    expect_equal(actual, expected)
  }
)

