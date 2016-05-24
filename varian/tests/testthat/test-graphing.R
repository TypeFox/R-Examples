context("Test Graphic Diagnostics")

set.seed(1234) # make reproducible
x <- matrix(rnorm(1000), ncol = 2)

test_that("vmp_plot returns a ggplot2 graph", {
  expect_that(vmp_plot(x, plot=FALSE)$Individual[[1]], is_a("ggplot"))
})
