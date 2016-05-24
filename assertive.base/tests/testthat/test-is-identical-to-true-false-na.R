test_that(
  "test is_identical_to_false with false returns true", 
  {
    expect_true(is_identical_to_false(FALSE))
  }
)

test_that(
  "test is_identical_to_false with a false vector returns false", 
  {
    expect_false(is_identical_to_false(logical(2)))
  }
)

test_that(
  "test is_identical_to_false with NA returns false", 
  {
    expect_false(is_identical_to_false(NA))
  }
)

test_that(
  "test is_identical_to_false with false + attr returns allow_attributes", 
  {
    x <- c(truth = FALSE)
    expect_false(is_identical_to_false(x))
    expect_true(is_identical_to_false(x, allow_attributes = TRUE))
  }
)


test_that(
  "test is_identical_to_na with TRUE returns false", 
  {
    expect_false(is_identical_to_na(TRUE))
  }
)

test_that(
  "test is_identical_to_na with NA returns true", 
  {
    expect_true(is_identical_to_na(NA))
  }
)

test_that(
  "test is_identical_to_na with an NA vector.returns false", 
  {
    expect_false(is_identical_to_na(rep.int(NA, 2)))
  }
)

test_that(
  "test is_identical_to_na NA + attr returns allow_attributes", 
  {
    x <- c(truth = NA)
    expect_false(is_identical_to_na(x))
    expect_true(is_identical_to_na(x, allow_attributes = TRUE))
  }
)


test_that(
  "test is_identical_to_true NA returns false", 
  {
    expect_false(is_identical_to_true(NA))
  }
)

test_that(
  "test is_identical_to_true with true returns true", 
  {
    expect_true(is_identical_to_true(TRUE))
  }
)

test_that(
  "test is_identical_to_true with a true vector returns false", 
  {
    expect_false(is_identical_to_true(rep.int(TRUE, 2)))
  }
)

test_that(
  "test is_identical_to_true true + attr returns allow_attributes", 
  {
    x <- c(truth = TRUE)
    expect_false(is_identical_to_true(x))
    expect_true(is_identical_to_true(x, allow_attributes = TRUE))
  }
)
