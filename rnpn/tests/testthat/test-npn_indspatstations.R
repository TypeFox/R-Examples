context("npn_indspatstations")

test_that("npn_indspatstations works well", {
  skip_on_cran()

  aa <- npn_indspatstations(speciesid = 35, stationid = c(60, 259), year = 2009)

  expect_is(aa, "data.frame")
  expect_is(aa$individual_id, "character")

  expect_named(aa, c('individual_id', 'individual_name', 'number_observations'))
  expect_gt(NROW(aa), 0)
})

test_that("when no match, returns NULL", {
  skip_on_cran()

  expect_null(npn_indspatstations(speciesid = "asdfasfasf", stationid = "Adfadsdf"))
})
