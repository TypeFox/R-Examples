context("npn_stationswithspp")

test_that("npn_stationswithspp works well", {
  skip_on_cran()

  aa <- npn_stationswithspp(speciesid = c(52,53,54))
  bb <- npn_stationswithspp(speciesid = 53)

  expect_is(aa, "data.frame")
  expect_is(aa$latitude, "character")
  expect_is(aa$station_name, "character")
  expect_gt(NROW(aa), 0)

  expect_is(bb, "data.frame")
  expect_is(bb$latitude, "character")
  expect_is(bb$station_name, "character")
  expect_gt(NROW(bb), 0)
})

test_that("when no match, returns NULL", {
  skip_on_cran()

  aa <- npn_stationswithspp(speciesid = "Adfadsdf")

  expect_null(aa)
})
