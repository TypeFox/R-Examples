context("map_gist works correctly")

test_that("map_gist works as expected", {
  skip_on_cran()

  library("spocc")

  spp <- c('Danaus plexippus','Accipiter striatus','Pinus contorta')
  dat <- occ(spp, from = 'gbif', limit = 30, has_coords = TRUE)
  dat <- fixnames(dat, "query")

  # Define colors
  g <- suppressMessages(map_gist(dat, color = c('#976AAE','#6B944D','#BD5945'), browse = FALSE))

  expect_is(g, "gist")
  expect_is(g$url, "character")
  expect_equal(g$description, "")
})
