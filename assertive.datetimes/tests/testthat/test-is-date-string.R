test_that(
  "test.is_date_string.a_character_vector.returns_true_when_string_contains_a_date", 
  {
    x <- c("1999-12-31 23:59:59", "1979-08-01 01:00:00", "31 Dec 1999 11:59:59PM", 
           "not a date", "NA")
    expected <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
    expect_equal(strip_attributes(actual <- is_date_string(x)), expected)
    expect_equal(names(actual), unname(x))
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "bad format"), c(2, 3)))
    )
  }
)
