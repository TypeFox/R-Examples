test_that("test.is_in_closed_range.0_to_4_in_1_to_3.returns_true_inside_bounds", 
{
  x <- c(0:4, NA)
  expected <- c(FALSE, TRUE, TRUE, TRUE, FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_in_closed_range(x, 1, 3)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("too low", "", "too high", "missing"), c(1, 3, 1, 1)))
  )
})

test_that("test.is_in_left_open_range.0_to_4_in_1_to_3.returns_true_inside_bounds", 
{
  x <- c(0:4, NA)
  expected <- c(FALSE, FALSE, TRUE, TRUE, FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_in_left_open_range(x, 1, 3)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("too low", "", "too high", "missing"), c(2, 2, 1, 1)))
  )
})

test_that("test.is_in_open_range.0_to_4_in_1_to_3.returns_true_inside_bounds", 
{
  x <- c(0:4, NA)
  expected <- c(FALSE, FALSE, TRUE, FALSE, FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_in_open_range(x, 1, 3)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("too low", "", "too high", "missing"), c(2, 1, 2, 1)))
  )
})

test_that("test.is_in_range.0_to_4_in_1_to_3.returns_true_inside_bounds", 
{
  x <- c(0:4, NA)
  expected <- c(FALSE, TRUE, TRUE, TRUE, FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_in_range(x, 1, 3)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("too low", "", "too high", "missing"), c(1, 3, 1, 1)))
  )
})

test_that("test.is_in_right_open_range.0_to_4_in_1_to_3.returns_true_inside_bounds", 
  {
    x <- c(0:4, NA)
    expected <- c(FALSE, TRUE, TRUE, FALSE, FALSE, NA)
    expect_equal(
      strip_attributes(actual <- is_in_right_open_range(x, 1, 3)), 
      expected
    )
    expect_named(actual)
    expect_equal(
      cause(actual),
      noquote(rep.int(c("too low", "", "too high", "missing"), c(1, 2, 2, 1)))
    )
  })

test_that("test.is_negative.minus_2_to_2.returns_true_when_negative", 
{
  x <- c(-2:2, NA)
  expected <- rep.int(c(TRUE, FALSE, NA), c(2, 3, 1))
  expect_equal(
    strip_attributes(actual <- is_negative(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "too high", "missing"), c(2, 3, 1)))
  )
})

test_that("test.is_non_negative.minus_2_to_2.returns_true_when_non_negative", 
{
  x <- c(2:-2, NA)
  expected <- rep.int(c(TRUE, FALSE, NA), c(3, 2, 1))
  expect_equal(
    strip_attributes(actual <- is_non_negative(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "too low", "missing"), c(3, 2, 1)))
  )
})

test_that("test.is_non_positive.minus_2_to_2.returns_true_when_non_positive", 
{
  x <- c(-2:2, NA)
  expected <- rep.int(c(TRUE, FALSE, NA), c(3, 2, 1))
  expect_equal(
    strip_attributes(actual <- is_non_positive(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "too high", "missing"), c(3, 2, 1)))
  )
})

test_that("test.is_percentage.minus_minus_1_to_101.returns_true_when_percentage", 
{
  x <- c(-1:101, NA)
  expected <- c(FALSE, rep.int(TRUE, 101), FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_percentage(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("too low", "", "too high", "missing"), c(1, 101, 1, 1)))
  )
})

test_that("test.is_positive.minus_2_to_2.returns_true_when_positive", 
{
  x <- c(2:-2, NA)
  expected <- rep.int(c(TRUE, FALSE, NA), c(2, 3, 1))
  expect_equal(
    strip_attributes(actual <- is_positive(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "too low", "missing"), c(2, 3, 1)))
  )
})

test_that("test.is_proportion.minus_minus_point_01_to_1_point_01.returns_true_when_percentage", 
{
  x <- c(seq.int(-0.01, 1.01, 0.01), NA)
  expected <- c(FALSE, rep.int(TRUE, 101), FALSE, NA)
  expect_equal(
    strip_attributes(actual <- is_proportion(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("too low", "", "too high", "missing"), c(1, 101, 1, 1)))
  )
}) 
