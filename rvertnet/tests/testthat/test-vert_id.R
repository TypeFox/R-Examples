context("vert_id")

test_that("vert_id works", {
  skip_on_cran()
  
  aa <- vert_id(ids = "urn:catalog:CM:Herps:116520", verbose = FALSE)
  
  expect_is(aa, "list")
  expect_is(aa$meta, "list")
  expect_is(aa$data, "data.frame")
  expect_equal(NROW(aa$data), 1)
  expect_named(aa$meta, c('request_date','response_records','request_origin','last_cursor',
                          'limit','query_version','matching_records','api_version'))
  expect_true(grepl("Bufo debilis", aa$data$scientificname))
})

test_that("vert_id works", {
  skip_on_cran()
  
  ids <- c("http://arctos.database.museum/guid/MSB:Mamm:56979?seid=1643089",
           "urn:catalog:CM:Herps:116520")
  aa <- vert_id(ids, verbose = FALSE)
  
  expect_is(aa, "list")
  expect_is(aa$meta, "list")
  expect_is(aa$data, "data.frame")
  expect_equal(NROW(aa$data), 2)
  expect_named(aa$meta, c('request_date','response_records','request_origin','last_cursor',
                          'limit','query_version','matching_records','api_version'))
  expect_true(any(grepl("Zapus", aa$data$scientificname)))
  expect_true(any(grepl("Bufo", aa$data$scientificname)))
})
