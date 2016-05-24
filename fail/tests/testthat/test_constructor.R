context("Constructor")

test_that("Creation of directories and checking of existing files", {
  path = tempfile()
  fail(path) # works on new dirs
  fail(path) # works on existing dirs
  path = tempfile()
  file.create(path)
  expect_error(fail(path))
})

test_that("Constructor checks input", {
  path = tempfile()
  expect_error(fail(path, extension="^"))
  expect_error(fail(path, extension=".RData"))
})

test_that("Argument 'all.files' works", {
  path = tempfile()
  f = fail(path, all.files = FALSE)
  expect_error(f$put(.x = 1), "hidden")

  f = fail(path, all.files = TRUE)
  f$put(.x = 1)
  expect_equal(f$ls(), ".x")
})
