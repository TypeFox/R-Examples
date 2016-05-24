test_that(
  "test.is_empty_character.a_character_vector.returns_true_when_string_is_missing_or_empty", 
  {
    x <- c(missing = NA_character_, empty = "", non_empty = "a", space = " ", 
      not_missing1 = "NA", not_missing2 = "<NA>")
    expected <- c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
    expect_equal(
      strip_attributes(actual <- is_empty_character(x)), 
      expected
    )
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("missing", "", "nonempty"), c(1, 1, 4)))
    )
  }
)

test_that(
  "test.is_non_empty_character.a_character_vector.returns_true_when_string_is_missing_or_empty", 
  {
    x <- c(missing = NA_character_, empty = "", non_empty = "a", space = " ", 
      not_missing1 = "NA", not_missing2 = "<NA>")
    expected <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
    expect_equal(
      strip_attributes(actual <- is_non_empty_character(x)), 
      expected
    )
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "empty", ""), c(1, 1, 4)))
    )
  }
)

test_that(
  "test.is_missing_or_empty_character.a_character_vector.returns_true_when_string_is_missing_or_empty", 
  {
    x <- c(missing = NA_character_, empty = "", non_empty = "a", space = " ", 
      not_missing1 = "NA", not_missing2 = "<NA>")
    expected <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    expect_equal(
      strip_attributes(actual <- is_missing_or_empty_character(x)), 
      expected
    )
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "nonempty"), c(2, 4)))
    )
  }
)

test_that(
  "test.is_non_missing_nor_empty_character.a_character_vector.returns_true_when_string_is_not_missing_nor_empty", 
  {
    x <- c(missing = NA_character_, empty = "", non_empty = "a", space = " ", 
      not_missing1 = "NA", not_missing2 = "<NA>")
    expected <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
    expect_equal(
      strip_attributes(actual <- is_non_missing_nor_empty_character(x)), 
      expected
    )
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("missing", "empty", ""), c(1, 1, 4)))
    )
  }
)
