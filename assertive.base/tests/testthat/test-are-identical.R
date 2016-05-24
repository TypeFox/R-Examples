test_that(
  "test are_identical with identical objects returns true", 
  {
    expect_true(are_identical((1:3) ^ 3, c(1, 8, 27)))
  }
)

test_that(
  "test are_identical with non-identical objects returns false", 
  {
    expect_false(are_identical((1:3) ^ 3, c(1, 8, 27.0000001)))
  }
)
