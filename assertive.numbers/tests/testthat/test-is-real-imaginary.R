test_that("test.is_imaginary.imaginary_numbers.returns_true_when_purely_imaginary", 
  {
    x <- c(0 + 0i, 0 + 1i, 1 + 0i, 1 + 1i, Inf, NA_complex_)
    expected <- rep.int(c(TRUE, FALSE, NA), c(2, 3, 1))
    expect_equal(
      strip_attributes(actual <- is_imaginary(x)), 
      expected
    )
    expect_named(actual)
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "real", "missing"), c(2, 3, 1)))
    )
  })

test_that("test.is_imaginary.real_numbers.returns_true_when_0", 
{
  x <- c(0, 1, -1, Inf, NA_real_)
  expected <- rep.int(c(TRUE, FALSE, NA), c(1, 3, 1))
  expect_equal(
    strip_attributes(actual <- is_imaginary(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "real", "missing"), c(1, 3, 1)))
  )
})

test_that("test.is_real.imaginary_numbers.returns_true_when_purely_real", 
{
  x <- c(0 + 0i, 1 + 0i, 0 + 1i, 1 + 1i, Inf * 1i, NA_complex_)
  expected <- rep.int(c(TRUE, FALSE, NA), c(2, 3, 1))
  expect_equal(
    strip_attributes(actual <- is_real(x)), 
    expected
  )
  expect_named(actual)
  expect_equal(
    cause(actual),
    noquote(rep.int(c("", "imaginary", "missing"), c(2, 3, 1)))
  )
})

test_that("test.is_real.real_numbers.returns_true_always", 
{
  x <- c(0, 1, -1, Inf, NA_real_)
  expected <- rep.int(TRUE, 5)
  expect_equal(
    strip_attributes(actual <- is_real(x)), 
    expected
  )
  expect_named(actual)
}) 
