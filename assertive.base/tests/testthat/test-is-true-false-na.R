test_that(
  "test is_false with a logical input returns true when false",
  {
    x <- c(TRUE, FALSE, NA)
    expected <- c(FALSE, TRUE, FALSE)
    actual <- assertive.base::is_false(x)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(cause(actual), noquote(c("true", "", "missing")))
  }
)

test_that(
  "test is_na with a logical input returns true when NA", 
  {
    x <- c(TRUE, FALSE, NA)
    expected <- c(FALSE, FALSE, TRUE)
    actual <- is_na(x)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(cause(actual), noquote(c("true", "false", "")))
  }
) 

test_that(
  "test is_na with a character input and no coercion returns true when not NA", 
  {
    x <- c("T", "F", "0", "1", "a", "NA", NA)
    expected <- rep.int(c(FALSE, TRUE), c(6, 1))
    actual <- is_na(x, coerce_to_logical = FALSE)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(
      cause(actual), 
      noquote(rep.int(c("not missing", ""), c(6, 1)))
    )
  }
) 

test_that(
  "test is_na with a character input and coercion returns true when not 'T' or 'F'", 
  {
    x <- c("T", "F", "0", "1", "a", "NA", NA)
    expected <- rep.int(c(FALSE, TRUE), c(2, 5))
    expect_warning(
      actual <- is_na(x, coerce_to_logical = TRUE),
      "Coercing x to class .logical.\\."
    )
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(
      cause(actual), 
      noquote(rep.int(c("true", "false", ""), c(1, 1, 5)))
    )
  }
) 

test_that(
  "test is_na with a complex input and no coercion returns true when not NA", 
  {
    x <- c(0, 1, 1i, 1 + 1i, NA)
    expected <- rep.int(c(FALSE, TRUE), c(4, 1))
    actual <- is_na(x, coerce_to_logical = FALSE)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(
      cause(actual), 
      noquote(rep.int(c("not missing", ""), c(4, 1)))
    )
  }
) 

test_that(
  "test is_na with a complex input and coercion returns true when not NA", 
  {
    x <- c(0, 1, 1i, 1 + 1i, NA)
    expected <- rep.int(c(FALSE, TRUE), c(4, 1))
    expect_warning(
      actual <- is_na(x, coerce_to_logical = TRUE),
      "Coercing x to class .logical.\\."
    )
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(
      cause(actual), 
      noquote(rep.int(c("false", "true", ""), c(1, 3, 1)))
    )
  }
) 

test_that(
  "test is_na with a list input and no coercion returns true when elements contain a single NA", 
  {
    x <- list(NA, TRUE, FALSE, c(TRUE, FALSE, NA))
    expected <- rep.int(c(TRUE, FALSE), c(1, 3))
    actual <- is_na(x, coerce_to_logical = FALSE)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(
      cause(actual), 
      noquote(rep.int(c("", "not missing"), c(1, 3)))
    )
  }
) 

test_that(
  "test is_na with a list input and coercion throw an error", 
  {
    x <- list(NA, TRUE, FALSE, c(TRUE, FALSE, NA))
    expect_error(
      is_na(x, coerce_to_logical = TRUE),
      "x cannot be coerced to any of these types: .logical.\\."
    )
  }
) 

test_that(
  "test is_not_false with a logical input returns true when not false", 
  {
    x <- c(TRUE, FALSE, NA)
    expected <- c(TRUE, FALSE, TRUE)
    actual <- is_not_false(x)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(cause(actual), noquote(c("", "false", "")))
  }
)

test_that(
  "test is_not_na with a logical input returns true when not NA", 
  {
    x <- c(TRUE, FALSE, NA)
    expected <- c(TRUE, TRUE, FALSE)
    actual <- is_not_na(x)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(cause(actual), noquote(c("", "", "missing")))
  }
) 

test_that(
  "test is_not_true with a logical input returns true when not true", 
  {
    x <- c(TRUE, FALSE, NA)
    expected <- c(FALSE, TRUE, TRUE)
    actual <- is_not_true(x)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(cause(actual), noquote(c("true", "", "")))
  }
) 

test_that(
  "test is_true with a logical input returns true when true", 
  {
    x <- c(TRUE, FALSE, NA)
    expected <- c(TRUE, FALSE, FALSE)
    actual <- assertive.base::is_true(x)
    expect_equal(strip_attributes(actual), expected)
    expect_named(actual)
    expect_equal(cause(actual), noquote(c("", "false", "missing")))
  }
) 
