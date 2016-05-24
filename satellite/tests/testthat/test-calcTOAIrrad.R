# devtools::test(".", "calcTOAIrrad")
context("Solar irradiation (ESun)")

#-------------------------------------------------------------------------------
test_that("calcTOAIrradModel for data frames works as expected", {
  calcTOAIrradModel(lut$L8_RSR, model = "MNewKur")
  calcTOAIrradModel(lut$L8_RSR, model = "MNewKur", normalize = FALSE, 
                    esd = calcEarthSunDist("2015-01-01"))
  calcTOAIrradModel(lut$L8_RSR, model = "MNewKur", normalize = FALSE, 
                    esd = calcEarthSunDist("2015-07-07"))
  calcTOAIrradModel(lut$L7_RSR, model = "MNewKur")
  calcTOAIrradModel(lut$L7_RSR, model = "MNewKur", normalize = FALSE, 
                    esd = calcEarthSunDist("2015-01-01"))
  calcTOAIrradModel(lut$L7_RSR, model = "MNewKur", normalize = FALSE, 
                    esd = calcEarthSunDist("2015-07-07"))
})


#-------------------------------------------------------------------------------
test_that("calcTOAIrradModel for satellite works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)  
  test <- calcTOAIrradModel(sat)
  
  expect_equal(as.character(getSatBID(test)[1]), "1")
  expect_equal(round(as.numeric(getSatESUN(test)[1]),4), round(1888.4115033, 4))
  expect_equal(as.character(getSatBID(test)[2]), "2")
  expect_equal(round(as.numeric(getSatESUN(test)[2]),4), round(1974.8429354, 4))
  expect_equal(as.character(getSatBID(test)[3]), "3")
  expect_equal(round(as.numeric(getSatESUN(test)[3]),4), round(1851.7520559, 4))
  expect_equal(as.character(getSatBID(test)[11]), "11")
  expect_equal(round(as.numeric(getSatESUN(test)[11]),4), round(0.1068904, 4))
  
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files[c(1,3,4)])  
  test <- calcTOAIrradModel(sat)
  
  expect_equal(as.character(getSatBID(test)[1]), "1")
  expect_equal(round(as.numeric(getSatESUN(test)[1]),4), round(1888.4115033, 4))
  expect_equal(as.character(getSatBID(test)[2]), "2")
  expect_equal(round(as.numeric(getSatESUN(test)[2]),4), round(1974.8429354, 4))
  expect_equal(as.character(getSatBID(test)[3]), "11")
  expect_equal(round(as.numeric(getSatESUN(test)[3]),4), round(0.1068904, 4))
})


#-------------------------------------------------------------------------------
test_that("calcTOAIrradRadRef for numeric works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)  
  test <- calcTOAIrradRadRef(x = getSatRadMax(sat, getSatBCDESolar(sat)), 
                             ref_max = getSatRefMax(sat, getSatBCDESolar(sat)), 
                             normalize = FALSE,
                             esd = calcEarthSunDist("2013-07-07")) 
  
  expect_equal(round(as.numeric(test[1]), 3), round( 1907.999, 3))
  expect_equal(round(as.numeric(test[3]), 3), round( 1800.423, 3))
})


#-------------------------------------------------------------------------------
test_that("calcTOAIrradRadRef for Satellite works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)  
  test <- calcTOAIrradRadRef(sat)
  
  expect_equal(as.character(getSatBID(test)[1]), "1")
  expect_equal(round(as.numeric(getSatESUN(test)[1]), 3), round(  1972.254, 3))
  expect_equal(as.character(getSatBID(test)[3]), "3")
  expect_equal(round(as.numeric(getSatESUN(test)[3]), 3), round( 1861.055, 3))
  
  sat <- satellite(files[c(1,3,4)])
  test <- calcTOAIrradRadRef(sat)
  
  expect_equal(as.character(getSatBID(test)[1]), "1")
  expect_equal(round(as.numeric(getSatESUN(test)[1]), 3), round( 1972.254, 3))
  expect_equal(as.character(getSatBID(test)[2]), "2")
  expect_equal(round(as.numeric(getSatESUN(test)[2]), 3), round( 2019.612, 3))
})


#-------------------------------------------------------------------------------
test_that("calcTOAIrradTable for character works as expected", {
  calcTOAIrradTable(x = "LE5")
  calcTOAIrradTable(x = "LE7")
  calcTOAIrradTable(x = "LE7", normalize = FALSE, 
                    esd = calcEarthSunDist("2015-01-01"))
})


#-------------------------------------------------------------------------------
test_that("calcTOAIrradTable for Satellite works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LE7*.tif"), full.names = TRUE)
  sat <- satellite(files)
  test <- calcTOAIrradTable(sat)
  
  expect_equal(as.character(getSatBID(test)[2]), "2")
  expect_equal(as.numeric(getSatESUN(test)[2]), 1812.0)
  
  sat <- satellite(files[c(1, 3, 6)])
  test <- calcTOAIrradTable(sat)
  
  expect_equal(as.character(getSatBID(test)[2]), "3")
  expect_equal(as.numeric(getSatESUN(test)[2]), 1533.0)
  
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)  
  expect_error(calcTOAIrradTable(sat), "Satellite ID LC8 is not supported, yet.")
})


#-------------------------------------------------------------------------------
test_that("Depricated satTOAIrrad for Satellite works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)  
  satTOAIrrad(sat, method = "Model")
  satTOAIrrad(sat, method = "RadRef")
})