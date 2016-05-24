context("Testing beta-utils \n")

test_that("beta names are correct", {
  expect_identical(get_beta_names("normal"), c("mu", "sigma"))
  expect_identical(get_beta_names("exp"), "lambda")
})