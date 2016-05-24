test_that("check vscat",{
  vscat(Petal.Length, Petal.Width, iris, "scat01", "iris", color=Species)
  expect_equal(file.exists("scat01.iris.svg.html"), TRUE)
})

test_that("check vhist",{
  vhist(Sepal.Length, iris, "hist01", "iris", fill=Species)
  expect_equal(file.exists("hist01.iris.svg.html"), TRUE)
})

test_that("check vpcp", {
  vpcp(iris, 1:4, "pcp01", "iris", groupColumn="Species")
  expect_equal(file.exists("pcp01.iris.svg.html"), TRUE)
})

test_that("check vcmap", {
  data(vsfuk2012)
  shp.path <- file.path(system.file(package="vdmR"), "etc/shapes/fukuoka2012.shp")
  frcol <- ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=median(vsfuk2012$FertilityRate))
  vcmap(shp.path, vsfuk2012, "CityCode", "CityCode", "map01", "vsfuk2012", fill=FertilityRate, ggscale=frcol)
  expect_equal(file.exists("map01.vsfuk2012.svg.html"), TRUE)
})

test_that("check vlaunch", {
  vlaunch(iris, "main", "iris", browse=FALSE)
  expect_equal(file.exists("main.iris.html"), TRUE)
})
