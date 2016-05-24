

context("examples for covr")
test_that("all basic examples", {

  a <- ggplot(dd_example("barchart"))
  a <- ggplot(dd_example("histogram"))
  a <- ggplot(dd_example("parcoord"))
  a <- ggplot(dd_example("plot1d"))
  a <- ggplot(dd_example("scatmat"))
  a <- ggplot(dd_example("timeseries"))
  a <- ggplot(dd_example("tour1d"))
  a <- ggplot(dd_example("tour2d"))
  a <- ggplot(dd_example("tour2x1d"))
  a <- ggplot(dd_example("xyplot"))

  expect_true(TRUE)
})
