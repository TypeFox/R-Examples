context("npn_stationsbystate")

test_that("npn_stationsbystate works well", {
  skip_on_cran()

  aa <- npn_stationsbystate()

  expect_is(aa, "data.frame")
  expect_named(aa, c('state', 'number_stations'))
  expect_is(aa$state, "character")
  expect_is(aa$number_stations, "character")
})
