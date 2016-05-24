context("directory functions")

test_that("isDirectory", {
  expect_true(isDirectory("."))
  expect_identical(isDirectory(".", ".."), c(TRUE, TRUE))
  expect_false(isDirectory("foofoo"))
  expect_identical(isDirectory(".", "foofoo"), c(TRUE, FALSE))
})


test_that("isEmptyDirectory", {
  expect_false(isEmptyDirectory("."))
  td = tempfile()
  dir.create(td)
  expect_true(isEmptyDirectory(td))
  expect_identical(isEmptyDirectory(td, ".."), c(TRUE, FALSE))
  expect_false(isEmptyDirectory("foofoo"))
  expect_false(isEmptyDirectory(tempfile()))
})



