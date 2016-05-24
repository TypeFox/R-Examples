context("map_ggmap works correctly")

test_that("map_ggmap maps - spocc input", {
  skip_on_cran()

  library("spocc")
  gd <- occ(query = 'Accipiter striatus', from = 'gbif', limit=25, has_coords = TRUE)
  ff <- suppressMessages(map_ggmap(gd))
  gg <- suppressMessages(map_ggmap(gd$gbif))

  expect_is(ff, "gg")
  expect_is(ff$coordinates, "CoordMap")
  expect_is(gg, "gg")
  expect_is(gg$coordinates, "CoordMap")
  expect_identical(ff$data, gg$data)
})


test_that("map_ggmap maps - rgbif input", {
  skip_on_cran()

  library("rgbif")
  res <- occ_search(scientificName = "Puma concolor", limit = 20)
  gg <- suppressMessages(map_ggmap(res))

  expect_is(gg, "gg")
  expect_is(gg$data, "data.frame")
  expect_true("lon" %in% names(gg$data))
  expect_match(gg$labels$title, "Puma concolor")
})


test_that("map_ggmap maps - dat.frame input", {
  skip_on_cran()

  df <- data.frame(name = c('Poa annua', 'Puma concolor', 'Foo bar'),
                   longitude = c(-120, -121, -123),
                   latitude = c(41, 42, 45), stringsAsFactors = FALSE)
  gg <- suppressMessages(map_ggmap(df))

  expect_is(gg, "gg")
  expect_is(gg$data, "data.frame")
  expect_true("lon" %in% names(gg$data))
})


test_that("map_ggmap maps - fails well", {
  skip_on_cran()
  expect_error(map_ggmap(), "\"x\" is missing")
  expect_error(map_ggmap(5), "does not support input of class")
  expect_error(map_ggmap("thig"), "does not support input of class")
  expect_error(map_ggmap(list()), "does not support input of class")
  expect_error(map_ggmap(matrix()), "does not support input of class")
  expect_error(map_ggmap(new.env()), "does not support input of class")
})
