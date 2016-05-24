test_that(
  "test are_set_equal with equal sets returns true",  
  {
    x <- 1:5
    y <- c(1, 3, 5, 4, 2)
    expect_true(are_set_equal(x, y))
  }
)

test_that(
  "test are_set_equal with different size sets returns false",  
  {
    x <- 1:5
    y <- c(1:4, 4)
    expect_false(actual <- are_set_equal(x, y))
    expect_equal(
      cause(actual),
      noquote("1:5 and c(1, 2, 3, 4) have different numbers of elements (5 versus 4).")
    )
  }
)

test_that(
  "test are_set_equal with unequal sets returns false",  
  {
    x <- 1:5
    y <- c(99, 3, 5, 4, 2)
    expect_false(actual <- are_set_equal(x, y))
    # R is doing something odd with references to variable names, so the
    # cause isn't "The element '1' in x is not in y."
    # Using force() doesn't seem to make a difference.
    expect_equal(
      cause(actual),
      noquote("The element '1' in 1:5 is not in c(99, 3, 5, 4, 2).")
    )
  }
)

test_that(
  "test is_subset with equal sets returns true",  
  {
    x <- 1:5
    expect_true(is_subset(x, x))
  }
)

test_that(
  "test is_subset with equal sets and strictly = TRUE returns false",  
  {
    x <- 1:5
    expect_false(is_subset(x, x, strictly = TRUE))
  }
)

test_that(
  "test is_subset with a subset returns true",  
  {
    x <- 1:5
    y <- 6:1
    expect_true(is_subset(x, y))
  }
)

test_that(
  "test is_subset with a non-subset returns false",  
  {
    x <- 1:5
    y <- 4:1
    expect_false(actual <- is_subset(x, y))
    expect_equal(
      cause(actual),
      noquote("The element '5' in x is not in y.")
    )
  }
)

test_that(
  "test is_superset with equal sets returns true",  
  {
    x <- 1:5
    expect_true(is_superset(x, x))
  }
)

test_that(
  "test is_superset with equal sets and strictly = TRUE returns false",  
  {
    x <- 1:5
    expect_false(is_superset(x, x, strictly = TRUE))
  }
)

test_that(
  "test is_superset with a superset returns true",  
  {
    x <- 1:6
    y <- 5:1
    expect_true(is_superset(x, y))
  }
)

test_that(
  "test is_superset with a non-superset returns false",  
  {
    x <- 1:4
    y <- 5:1
    expect_false(actual <- is_superset(x, y))
    expect_equal(
      cause(actual),
      noquote("The element '5' in y is not in x.")
    )
  }
)

