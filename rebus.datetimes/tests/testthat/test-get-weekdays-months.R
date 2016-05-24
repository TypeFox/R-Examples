context("get_weekdays")

test_that(
  "get_weekdays returns a group of days of the week",
  {
    a_monday <- as.Date("2001-01-01")
    expected <- as.regex("(?:Monday|Tuesday|Wednesday|Thursday|Friday|Saturday|Sunday)")
    actual <- get_weekdays(from = a_monday)
    expect_equal(actual, expected)
  }
)

test_that(
  "get_weekdays with abbreviate = TRUE returns a group of days of the week",
  {
    a_monday <- as.Date("2001-01-01")
    expected <- as.regex("(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun)")
    actual <- get_weekdays(TRUE, from = a_monday)
    expect_equal(actual, expected)
  }
)
context("get_months")

test_that(
  "get_months returns a group of months of the year",
  {
    in_january <- as.Date("2001-01-01")
    expected <- as.regex("(?:January|February|March|April|May|June|July|August|September|October|November|December)")
    actual <- get_months(from = in_january)
    expect_equal(actual, expected)
  }
)

test_that(
  "get_months with abbreviate = TRUE returns a group of months of the year",
  {
    in_january <- as.Date("2001-01-01")
    expected <- as.regex("(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)")
    actual <- get_months(TRUE, from = in_january)
    expect_equal(actual, expected)
  }
)
