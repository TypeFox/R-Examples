library(ggplot2)
library(grid)

context("Plot: autoplot")
# Test autoplot.evalcurve
#

test_that("autoplot.evalcurve", {
  pdf(NULL)
  on.exit(dev.off())

  toolset <- create_toolset(set_names = "crv5")
  testset <- create_testset("curve", "c2")
  evalcrv1 <- run_evalcurve(testset, toolset)

  expect_error(suppressWarnings(autoplot(evalcrv1)), NA)

  toolset <- create_toolset(c("ROCR", "precrec"))
  testset <- create_testset("curve", "c2")
  evalcrv2 <- run_evalcurve(testset, toolset)

  expect_error(suppressWarnings(autoplot(evalcrv2)), NA)
})

test_that("autoplot.evalcurve ret_grob", {
  pdf(NULL)
  on.exit(dev.off())

  toolset <- create_toolset(set_names = "crv5")
  testset <- create_testset("curve", "c2")
  evalcrv <- run_evalcurve(testset, toolset)

  pp <- suppressWarnings(autoplot(evalcrv, ret_grob = TRUE))
  expect_true(is(pp, "grob"))

  expect_error(suppressWarnings(grid.draw(pp)), NA)

  expect_error(suppressWarnings(autoplot(evalcrv, ret_grob = 1)),
               "ret_grob is not a flag")
})

test_that("autoplot.evalcurve base_plot", {
  pdf(NULL)
  on.exit(dev.off())

  toolset <- create_toolset(set_names = "crv5")
  testset <- create_testset("curve", "c2")
  evalcrv <- run_evalcurve(testset, toolset)

  pp1 <- suppressWarnings(autoplot(evalcrv, base_plot = TRUE, ret_grob = TRUE))
  expect_equal(length(pp1$grobs), 6)

  pp2 <- suppressWarnings(autoplot(evalcrv, base_plot = FALSE, ret_grob = TRUE))
  expect_equal(length(pp2$grobs), 5)

  expect_error(suppressWarnings(autoplot(evalcrv, base_plot = 1)),
               "base_plot is not a flag")
})

test_that("autoplot.evalcurve ncol & nrow", {
  pdf(NULL)
  on.exit(dev.off())

  toolset <- create_toolset(set_names = "crv5")
  testset <- create_testset("curve", "c2")
  evalcrv <- run_evalcurve(testset, toolset)

  expect_error(suppressWarnings(autoplot(evalcrv, ncol = 1)),
               "Both ncol and nrow must be set")

  expect_error(suppressWarnings(autoplot(evalcrv, nrow = 1)),
               "Both ncol and nrow must be set")

  expect_error(suppressWarnings(autoplot(evalcrv, ncol = 1, nrow = 0)),
               "nrow not greater than 0")

  expect_error(suppressWarnings(autoplot(evalcrv, ncol = 0, nrow = 1)),
               "ncol not greater than 0")
})
