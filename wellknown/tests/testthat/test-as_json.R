context("as_json")

str <- "POLYGON ((100 0.1, 101.1 0.3, 101 0.5, 100 0.1),
   (103.2 0.2, 104.8 0.2, 100.8 0.8, 103.2 0.2))"

test_that("as_json works", {
  expect_is(wkt2geojson(str), "geojson")
  expect_is(unclass(wkt2geojson(str)), "list")

  expect_is(as_json(wkt2geojson(str)), "json")
})

test_that("as_json fails correctly", {
  expect_error(as_json(-116.4), "no applicable method")
  expect_error(as_json("a"), "no applicable method")
})
