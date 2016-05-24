context("npn_indsatstations")

test_that("npn_indsatstations works well", {
  skip_on_cran()

  aa <- npn_indsatstations(stationid = c(507, 523))

  expect_is(aa, "data.frame")
  expect_is(aa$individual_id, "character")
  expect_is(aa$kingdom, "character")
  expect_gt(NROW(aa), 0)

  expect_equal(unique(aa$kingdom), "Plantae")
})

test_that("when no match, returns NULL", {
  skip_on_cran()

  expect_null(npn_indsatstations(stationid = "Adfadsdf"))
  expect_null(npn_indsatstations(stationid = 1111111111111))
})
