context("map_ggplot works correctly")

test_that("map_ggplot work as expected", {
  skip_on_cran()

  library("ggplot2")
  library("spocc")

	res <- occ(query = "Lynx rufus californicus", from = "gbif")
  map1 <- suppressMessages(map_ggplot(res))
 	res2 <- occ(query = 'Accipiter striatus', from = 'gbif')
	map2 <- suppressMessages(map_ggplot(res2))
	expect_is(res, "occdat")
	expect_is(res2, "occdat")
	expect_is(map1, "ggplot")
	expect_is(map2, "ggplot")
	unlink("ggmapTemp.png")
})
