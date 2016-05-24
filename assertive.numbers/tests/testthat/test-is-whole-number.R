test_that(
  "test.is_whole_number.a_numeric_vector.returns_true_for_whole_numbers", 
  {
    x <- c(
      0, 1, 100 * .Machine$double.eps, 100 * -.Machine$double.eps, 
      -0.5, 101 * .Machine$double.eps, -101 * .Machine$double.eps, 
      Inf, -Inf, NaN, NA
    )
    expected <- rep.int(c(TRUE, FALSE, NA), c(4, 5, 2))
    expect_equal(
      strip_attributes(actual <- is_whole_number(x)), 
      expected
    )
    expect_named(actual)
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "fractional", "infinite", "missing"), c(4, 3, 2, 2)))
    )
  }
)

test_that(
  "test.is_whole_number.no_tolerance.returns_true_for_exactly_whole_numbers", 
  {
    x <- c(
      0, 1, 100 * .Machine$double.eps, 100 * -.Machine$double.eps, 
      -0.5, 101 * .Machine$double.eps, -101 * .Machine$double.eps, 
      Inf, -Inf, NaN, NA
    )
    expected <- rep.int(c(TRUE, FALSE, NA), c(2, 7, 2))
    expect_equal(
      strip_attributes(actual <- is_whole_number(x, 0)), 
      expected
    )
    expect_named(actual)
    expect_equal(
      cause(actual),
      noquote(rep.int(c("", "fractional", "infinite", "missing"), c(2, 5, 2, 2)))
    )
  }
) 
