context("lazy.citation")

test_that("markdown - lazy.citation",
{
  options(lazyReportFormat = "markdown")
  expect_equal(
    lazy.citation(),
    paste0("R Core Team, \"R: A Language and Environment for Statistical Computing,\" R Foundation for Statistical Computing, Vienna, Austria, Vol. , (", R.Version()$year, ") ."))
})

test_that("latex - lazy.citation",
{
  options(lazyReportFormat = "latex")
  expect_equal(
    lazy.citation(),
    paste0("R Core Team, ``R: A Language and Environment for Statistical Computing,\" R Foundation for Statistical Computing, Vienna, Austria, Vol. , (", R.Version()$year, ") ."))
})

test_that("latex - lazy.citation",
{
  options(lazyReportFormat = "html")
  expect_equal(
    lazy.citation(),
    paste0("R Core Team, \"R: A Language and Environment for Statistical Computing,\" R Foundation for Statistical Computing, Vienna, Austria, Vol. , (", R.Version()$year, ") ."))
})