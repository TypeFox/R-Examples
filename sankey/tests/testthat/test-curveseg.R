
context("Curveseg")

test_that("curve_seg is good", {

  ## This is a rather minimal test, but until we have some
  ## image fingerprinting on CRAN that works reliable, it
  ## is hard to do better

  mypolygon <- function() {
    num <- 0
    function(...) num <<- num + 1
  }

  mylines   <- function() {
    num <- 0
    function(...) num <<- num + 1
  }

  with_mock(
    `graphics::polygon` = mypolygon_func <- mypolygon(),
    `graphics::lines`   = mylines_func   <- mylines(),
    curveseg(0, 0, 1, 1, colorstyle = "color")
  )

  expect_equal(mypolygon_func(), 50)
  expect_equal(mylines_func(),   99)

  with_mock(
    `graphics::polygon` = mypolygon_func <- mypolygon(),
    `graphics::lines`   = mylines_func   <- mylines(),
    curveseg(0, 0, 1, 1, colorstyle = "color", curvestyle = "line")
  )

  expect_equal(mypolygon_func(), 50)
  expect_equal(mylines_func(),   99)
})
