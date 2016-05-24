context("map_ggplot works correctly")

test_that("map_ggplot work as expected", {
  skip_on_cran()

  library("ggplot2")
  library("spocc")

	ecoengine_data <- occ(query = "Lynx rufus californicus", from = "ecoengine")
  map1 <- suppressMessages(map_ggplot(ecoengine_data))
 	gbif_data <- occ(query = 'Accipiter striatus', from = 'gbif')
	map2 <- suppressMessages(map_ggplot(gbif_data))
	expect_is(ecoengine_data, "occdat")
	expect_is(gbif_data, "occdat")
	expect_is(map1, "ggplot")
	expect_is(map2, "ggplot")
	unlink("ggmapTemp.png")
})
