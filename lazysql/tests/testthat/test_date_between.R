context("Test date_between")

test_that("date_between works as expected", {
  date1 <- as.Date("2016-02-22")
  date2 <- as.Date("2016-02-11")
  date_column <- "STD_1"

  expect_identical(
    date_between(date_column, date1),
    "STD_1 between to_date('2016-02-22', 'yyyy-mm-dd') and to_date('2016-02-22', 'yyyy-mm-dd')")
  expect_identical(
    date_between(date_column, c(date1, date2)),
    "STD_1 between to_date('2016-02-11', 'yyyy-mm-dd') and to_date('2016-02-22', 'yyyy-mm-dd')")
})

test_that("date_between checks dates", {
  date1 <- as.Date("2016-02-22")
  date_column <- "STD"

  expect_error(date_between(date_column, as.POSIXct(date1)))
  expect_error(date_between(date_column, as.Date(NA)))
  expect_error(date_between(date_column, date1 + 1:11))
  expect_error(date_between(date_column, date1[0]))
  expect_error(date_between(date_column, NULL))
})

test_that("date_between checks column_names", {
  date1 <- as.Date("2016-02-22")

  expect_error(date_between("", date1))
  expect_error(date_between(NA, date1))
  expect_error(date_between(1L, date1))
  expect_error(date_between("'wrong'", date1))
  expect_error(date_between("wrong wrong", date1))
  expect_error(date_between("123", date1))
  expect_error(date_between("ABC$", date1))
})


