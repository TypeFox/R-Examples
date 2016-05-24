# TODO: tests for is_single_character

test_that(
  "test.is_numeric_string.a_character_vector.returns_true_when_string_contains_a_number", 
  {
    x <- c("1", "-2.3e4", "Inf", "one", "NA")
    expected <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
    expect_equal(
      strip_attributes(actual <- is_numeric_string(x)), 
      expected
    )
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "bad format"), c(3, 2)))
    )
  }
)


test_that(
  "test.is_logical_string.a character vector.returns true when string contains a logical value", 
  {
    x <- c(
      "TRUE", "FALSE", "true", "false", "True", "False", "T", "F", 
      "trUE", "FaLsE", "t", "f", "NA"
    )
    expected <- rep.int(c(TRUE, FALSE), c(8, 5))
    expect_equal(
      strip_attributes(actual <- is_logical_string(x)), 
      expected
    )
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "bad format"), c(8, 5)))
    )
  }
)

