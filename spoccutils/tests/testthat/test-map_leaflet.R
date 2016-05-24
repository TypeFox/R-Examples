context("map_leaflet works correctly")

test_that("Leaflet maps and geoJSON work", {
  skip_on_cran()

  library("spocc")
  spp <- c('Danaus plexippus','Accipiter striatus','Pinus contorta','Puma concolor', 'Ursus americanus','Gymnogyps californianus')
  dat <- occ(query = spp, from = 'gbif', limit = 25, gbifopts = list(hasCoordinate = TRUE))
  map_leaflet(dat, map_provider = 'toner', dest = ".")
  expect_true(file.exists("map"))
  expect_true(file.exists("data.geojson"))
  unlink("data.geojson")
  unlink("map/", recursive = TRUE)
})
