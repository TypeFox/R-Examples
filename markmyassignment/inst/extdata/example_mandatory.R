### Assignment structure ###

context("Mandatory tests")

test_that("Mandatory tests", {
  expect_true(exists("my_name"), "Variable my_name is missing")
  res <- use_package()("cheating_package")
  expect_that(res$passed, is_false(), info = "package 'cheating_package' should not be used.")
})

