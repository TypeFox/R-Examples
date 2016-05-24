test_that("test.is_comma_for_decimal_point.any_locale.returns_true_if_locale_uses_comma", 
{
  dp <- unname(Sys.localeconv()["decimal_point"])
  expected <- if(is.null(dp) || !nzchar(dp)) NA else dp == ","
  actual <- is_comma_for_decimal_point()
  expect_equal(strip_attributes(actual), expected)
  if (!isTRUE(actual)) 
  {
    msg <- if(is.null(dp))
    {
      "R has been compiled without support for locales."
    } else if(!nzchar(dp))
    {
      "The locale convention for a (numeric) decimal point has not been defined."
    } else
    {
      "The locale convention is to use a '.' for a (numeric) decimal point."
    }
    expect_equal(cause(actual), noquote(msg))
  }
})

test_that("test.is_comma_for_decimal_point.any_locale_money_type.returns_true_if_locale_uses_comma", 
{
  dp <- unname(Sys.localeconv()["mon_decimal_point"])
  expected <- if(is.null(dp) || !nzchar(dp)) NA else dp == ","
  actual <- is_comma_for_decimal_point("money")
  expect_equal(strip_attributes(actual), expected)
  if (!isTRUE(actual)) 
  {
    msg <- if(is.null(dp))
    {
      "R has been compiled without support for locales."
    } else if(!nzchar(dp))
    {
      "The locale convention for a (monetary) decimal point has not been defined."
    } else
    {
      "The locale convention is to use a '.' for a (monetary) decimal point."
    }
    expect_equal(cause(actual), noquote(msg))
  }
})

test_that("test.is_period_for_decimal_point.any_locale.returns_true_if_locale_uses_period", 
{
  dp <- unname(Sys.localeconv()["decimal_point"])
  expected <- if(is.null(dp) || !nzchar(dp)) NA else dp == "."
  actual <- is_period_for_decimal_point()
  expect_equal(strip_attributes(actual), expected)
  if (!isTRUE(actual)) 
  {
    msg <- if(is.null(dp))
    {
      "R has been compiled without support for locales."
    } else if(!nzchar(dp))
    {
      "The locale convention for a (numeric) decimal point has not been defined."
    } else
    {
      "The locale convention is to use a ',' for a (numeric) decimal point."
    }
    expect_equal(cause(actual), noquote(msg))
  }
})

test_that("test.is_period_for_decimal_point.any_locale_money_type.returns_true_if_locale_uses_period", 
{
  dp <- unname(Sys.localeconv()["mon_decimal_point"])
  expected <- if(is.null(dp) || !nzchar(dp)) NA else dp == "."
  actual <- is_period_for_decimal_point("money")
  expect_equal(strip_attributes(actual), expected)
  if (!isTRUE(actual)) 
  {
    msg <- if(is.null(dp))
    {
      "R has been compiled without support for locales."
    } else if(!nzchar(dp))
    {
      "The locale convention for a (monetary) decimal point has not been defined."
    } else
    {
      "The locale convention is to use a ',' for a (monetary) decimal point."
    }
    expect_equal(cause(actual), noquote(msg))
  }
})

