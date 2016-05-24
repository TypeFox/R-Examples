test_that("result of function", {
  expect_equal(class(deputies_get_ids("sejmrp", "reader", "qux94874", "services.mini.pw.edu.pl", 8, .Platform$OS.type == "windows")),
    "character")
})

test_that("names", {
  expect_named(deputies_get_ids("sejmrp", "reader", "qux94874", "services.mini.pw.edu.pl", 7, .Platform$OS.type == "windows"))
})


test_that("length", {
  expect_more_than(length(deputies_get_ids("sejmrp", "reader", "qux94874", "services.mini.pw.edu.pl", 7, .Platform$OS.type == "windows")), 0)
})