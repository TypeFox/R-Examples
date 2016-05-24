context("map_leaflet works correctly")

test_that("Leaflet maps and geoJSON work", {
  skip_on_cran()

  library("spocc")
  spp <- c('Danaus plexippus','Accipiter striatus','Pinus contorta','Puma concolor', 'Ursus americanus','Gymnogyps californianus')
  dat <- occ(query = spp, from = 'gbif', limit = 25, gbifopts = list(hasCoordinate = TRUE))
  gg <- map_leaflet(dat)
  expect_is(gg, "leaflet")
  expect_is(attr(gg$x, "leafletData"), "data.frame")
  expect_true("longitude" %in% names(attr(gg$x, "leafletData")))
})
