context("Swedish tests")

test_that("toOrdinal correctly processes integers 1-40 in Swedish", {
  first_40 <- c("1:a", "2:a", "3:e", "4:e", "5:e", "6:e", "7:e", "8:e", "9:e",
    "10:e", "11:e", "12:e", "13:e", "14:e", "15:e", "16:e", "17:e",
    "18:e", "19:e", "20:e", "21:a", "22:a", "23:e", "24:e", "25:e",
    "26:e", "27:e", "28:e", "29:e", "30:e", "31:a", "32:a", "33:e",
    "34:e", "35:e", "36:e", "37:e", "38:e", "39:e", "40:e")
  using_toOrdinal <- sapply(c(1:40), "toOrdinal", language="swedish")

  expect_equal(
    first_40, using_toOrdinal
  )
})

test_that("toOrdinal correctly processes integers 101, 102, 111, and 112 in Swedish", {
  adtl_nums <- c("101:a", "102:a", "111:e", "112:e")
  using_toOrdinal <- sapply(c(101, 102, 111, 112), "toOrdinal", language="swedish")

  expect_equal(
    adtl_nums, using_toOrdinal
  )
})

test_that("toOrdinal correctly errors when given a negative integer.",{

  expect_error(
    toOrdinal(-1), "Number supplied to 'toOrdinal' must be a positive integer."
  )
})
