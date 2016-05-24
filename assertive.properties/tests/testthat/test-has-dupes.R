test_that(
  "test.has_no_duplicates.without_duplicates.returns_false",
  {
    x <- 1:10
    expect_true(has_no_duplicates(x))
  }
)

test_that(
  "test.has_no_duplicates.with_1_duplicate.returns_false",
  {
    x <- rep.int(1, 2)
    actual <- has_no_duplicates(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("x has a duplicate at position 2."))
  }
)

test_that(
  "test.has_no_duplicates.with_multiple_duplicates.returns_false",
  {
    x <- rep.int(1, 3)
    actual <- has_no_duplicates(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("x has duplicates at positions 2, 3."))
  }
)


test_that(
  "test.has_duplicates.without_duplicates.returns_false",
  {
    x <- 1:10
    actual <- has_duplicates(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("x has no duplicates."))
  }
)

test_that(
  "test.has_duplicates.with_duplicates.returns_false",
  {
    x <- rep.int(1, 2)
    expect_true(has_duplicates(x))
  }
)

