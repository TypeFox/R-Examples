context("npn_obsspbyday")

test_that("npn_obsspbyday works well", {
  skip_on_cran()

  aa <- npn_obsspbyday(speciesid=357, startdate='2010-04-01', enddate='2012-01-05')

  expect_is(aa, "list")
  expect_named(aa, "357")
  expect_is(aa$`357`, "data.frame")
  expect_is(aa$`357`$date, "character")
  expect_is(aa$`357`$count, "character")

  expect_named(aa$`357`, c('date', 'count'))
})

test_that("when no match, returns empty data.frame", {
  skip_on_cran()

  expect_equal(NROW(npn_obsspbyday("Adfadsdf")[[1]]), 0)
})
