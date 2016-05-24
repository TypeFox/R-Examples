test_that(
  "test.has_cols.with_columns.returns_true",
  {
    x <- matrix(1:12, nrow = 3)
    expect_true(has_cols(x))
  }
)

test_that(
  "test.has_cols.without_columns.returns_true",
  {
    x <- 1:10
    actual <- has_cols(x)
    expect_false(actual)
    expect_equal(
      cause(actual), 
      noquote("The number of columns in x is NULL.")
    )
  }
)

test_that(
  "test.has_cols.zero_cols.returns_false",
  {
    x <- matrix(numeric(), ncol = 0)
    actual <- has_cols(x)
    expect_false(actual)
    expect_equal(
      cause(actual), 
      noquote("The number of columns in x is zero.")
    )
  }
)


test_that(
  "test.has_dims.a_matrix.returns_true",
  {
    mat <- matrix(1:12, nrow = 3)
    expect_true(has_dims(mat))
  }
)

test_that(
  "test.has_dims.a_data_frame.returns_true",
  {
    dfr <- data.frame(x = 1:5, y = runif(5))
    expect_true(has_dims(dfr))
  }
)

test_that(
  "test.has_dims.a_vector.returns_false",
  {
    x <- 1:3
    actual <- has_dims(x)
    expect_false(actual)
    expect_equal(cause(actual), noquote("The dimensions of x are NULL."))
  }
)


test_that(
  "test.has_rows.with_rows.returns_true",
  {
    x <- matrix(1:12, nrow = 3)
    expect_true(has_rows(x))
  }
)

test_that(
  "test.has_rows.without_rows.returns_false",
  {
    x <- 1:10
    actual <- has_rows(x)
    expect_false(actual)
    expect_equal(
      cause(actual), 
      noquote("The number of rows in x is NULL.")
    )
  }
)

test_that(
  "test.has_rows.zero_rows.returns_false",
  {
    x <- matrix(numeric(), nrow = 0)
    actual <- has_rows(x)
    expect_false(actual)
    expect_equal(
      cause(actual), 
      noquote("The number of rows in x is zero.")
    )
  }
)
