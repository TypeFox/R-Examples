# devtools::test(".", "satInfo")
context("satInfo")

test_that("satAddLog works as expected for Landsat 8", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LC8*.tif"), 
                      full.names = TRUE)
  sat <- satellite(files)
  
  sat <- addSatLog(sat)
  sat <- addSatLog(sat)
  expect_equal(length(getSatLog(sat)), 3)
  expect_equal(names(getSatLog(sat)[3]), "ps0003")
  
  #   sort( sapply(mget(ls()),object.size) )
  #   print(object.size(x=lapply(ls(), get)), units="Mb")
  
})