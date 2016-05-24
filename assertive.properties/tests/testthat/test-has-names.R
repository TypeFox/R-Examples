test_that(
  "test.has_colnames.with_colnames.returns_true",
  {
    x <- matrix(1:12, nrow = 3, dimnames = list(letters[1:3], LETTERS[1:4]))
    expect_true(has_colnames(x))
  }
)

test_that(
  "test.has_colnames.without_colnames.returns_false",
  {
    x <- matrix(1:12, nrow = 3)
    actual <- has_colnames(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The column names of x are NULL."))
  }
)

test_that(
  "test.has_colnames.with_empty_colnames.returns_false",
  {
    x <- matrix(
      1:12, 
      nrow = 3, 
      dimnames = list(character(3), character(4))
    )
    actual <- has_colnames(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The column names of x are all empty."))
  }
)


test_that(
  "test.has_dimnames.with_dimnames.returns_true",
  {
    x <- matrix(1:12, nrow = 3, dimnames = list(letters[1:3], LETTERS[1:4]))
    expect_true(has_dimnames(x))
  }
)

test_that(
  "test.has_dimnames.without_dimnames.returns_false",
  {
    x <- matrix(1:12, nrow = 3)
    actual <- has_dimnames(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The dimension names of x are NULL."))
  }
)

test_that(
  "test.has_dimnames.with_empty_dimnames.returns_false",
  {
    x <- matrix(
      1:12, 
      nrow = 3, 
      dimnames = list(character(3), character(4))
    )
    actual <- has_dimnames(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The dimension names of x are all empty."))
  }
)


test_that(
  "test.has_names.named_vector.returns_true",
  {
    x <- c(foo = 1, 2, 3)
    expect_true(assertive.properties::has_names(x))
  }
)

test_that(
  "test.has_names.data_frame.returns_true",
  {
    dfr <- data.frame(x = 1:5, y = runif(5))
    expect_true(assertive.properties::has_names(dfr))
  }
)

test_that(
  "test.has_names.unnamed_vector.returns_false",
  {
    x <- 1:3
    actual <- assertive.properties::has_names(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The names of x are NULL."))
  }
)


test_that(
  "test.has_rownames.with_rownames.returns_true",
  {
    x <- matrix(
      1:12, 
      nrow = 3, 
      dimnames = list(letters[1:3], LETTERS[1:4])
    )
    expect_true(has_rownames(x))
  }
)

test_that(
  "test.has_rownames.without_rownames.returns_false",
  {
    x <- matrix(1:12, nrow = 3)
    actual <- has_rownames(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The row names of x are NULL."))
  }
)

test_that(
  "test.has_rownames.with_empty_rownames.returns_false",
  {
    x <- matrix(
      1:12, 
      nrow = 3, 
      dimnames = list(character(3), character(4))
    )
    actual <- has_rownames(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The row names of x are all empty."))
  }
)
