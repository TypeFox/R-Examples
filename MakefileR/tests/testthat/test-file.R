context("file")

test_that("Empty Makefile includes only comment", {
  makefile <- format(makefile())
  expect_match(makefile, "^# .* do not edit.*$")
  expect_equal(length(makefile), 1)
})

test_that("Printing works as expected", {
  with_mock(
    cat = function(x, sep) x,
    makefile <- print(makefile()))
  expect_match(makefile, "^# .* do not edit.*\n$")
  expect_equal(length(makefile), 1)
})
