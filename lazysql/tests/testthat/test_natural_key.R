context("Test natural_key")

test_that("natural_key works as expected", {
  expect_identical(
    natural_key(c("TAB1", "tab_2"),c("COL1", "col_2")),
    "TAB1.COL1 = tab_2.COL1 and TAB1.col_2 = tab_2.col_2")
})

test_that("natural_key checks table names", {
  expect_error(natural_key(  "TAB1",                "COL1"))
  expect_error(natural_key(c("TAB1", "TAB1"),       "COL1"))
  expect_error(natural_key(c("TAB1", ""),           "COL1"))
  expect_error(natural_key(c("TAB1", NA),           "COL1"))
  expect_error(natural_key(c("TAB1", 1L),           "COL1"))
  expect_error(natural_key(c("TAB1", "'noHyphen'"), "COL1"))
  expect_error(natural_key(c("TAB1", "no blanks"),  "COL1"))
  expect_error(natural_key(c("TAB1", "123"),        "COL1"))
  expect_error(natural_key(c("TAB1", "ABC$"),       "COL1"))
})

test_that("natural_key checks column_names", {
  expect_error(natural_key("TAB1", NULL))
  expect_error(natural_key("TAB1", character(0)))
  expect_error(natural_key("TAB1", character(1)))
  expect_error(natural_key("TAB1", ""))
  expect_error(natural_key("TAB1", NA))
  expect_error(natural_key("TAB1", 1L))
  expect_error(natural_key("TAB1", "'noHyphen'"))
  expect_error(natural_key("TAB1", "no blanks"))
  expect_error(natural_key("TAB1", "123"))
  expect_error(natural_key("TAB1", "ABC$"))
})


