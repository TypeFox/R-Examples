# devtools::test(".", "satInfo")
context("convDN2RU")


#-------------------------------------------------------------------------------
test_that("convDN2RU for Satellite works as expected", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LC8*.tif"), 
                      full.names = TRUE)
  sat <- satellite(files)
  
  t1 <- convDN2RU(sat, convert = "Rad")
  t2 <- convDN2RU(sat, convert = "Ref", szen_correction = TRUE)
  t3 <- convDN2RU(sat, convert = "Ref", szen_correction = FALSE)
  t4 <- convDN2RU(sat, convert = "BT")
  t5 <- convDN2RU(sat, convert = "all")
  
  expect_equal(
    round(raster::getValues(getSatDataLayer(t1, bcde = "B002n_RAD"))[50],4), 
    round(62.95932, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t2, bcde = "B002n_REF"))[50],4), 
    round( 0.1179185, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t3, bcde = "B002n_REF_NoSZEN"))[50],4), 
    round(0.10124, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t4, bcde = "B010n_BT"))[50],4), 
    round(303.5775, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t5, bcde = "B009n_RAD"))[50],4), 
    round(0.1295196, 4))
  
  
  bcde <- "B002n"
  
  t1 <- convDN2RU(x = getSatDataLayer(sat, bcde),
                        mult = getSatREFM(sat, bcde),
                        add = getSatREFA(sat, bcde))
  
  t2 <- convDN2RU(x = getSatDataLayer(sat, bcde),
                        mult = getSatREFM(sat, bcde),
                        add = getSatREFA(sat, bcde),
                        szen = getSatSZEN(sat, bcde))
  
  t3 <- convDN2RU(x = getSatDataLayer(sat, bcde),
                        mult = getSatRADM(sat, bcde),
                        add = getSatRADA(sat, bcde))
  
  bcde <- "B010n"
  t4 <- convDN2RU(x = getSatDataLayer(sat, bcde),
                        mult = getSatRADM(sat, bcde),
                        add = getSatRADA(sat, bcde),
                        k1 = getSatBTK1(sat, bcde),
                        k2 = getSatBTK2(sat, bcde))
  
  expect_equal(round(raster::getValues(t1)[50],3), round(0.10124,3))
  expect_equal(round(raster::getValues(t2)[50],3), round(0.1179185,3))
  expect_equal(round(raster::getValues(t3)[50],3), round(62.959,3))
  expect_equal(round(raster::getValues(t4)[50],3), round(303.5775,3))
  
})


#-------------------------------------------------------------------------------
test_that("Depricated satCalib and calibLinear works as expected", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LC8*.tif"), 
                      full.names = TRUE)
  sat <- satellite(files)
  
  t1 <- satCalib(sat, convert = "Rad")
  t2 <- satCalib(sat, convert = "Ref", szen_correction = TRUE)
  t3 <- satCalib(sat, convert = "Ref", szen_correction = FALSE)
  t4 <- satCalib(sat, convert = "BT")
  t5 <- satCalib(sat, convert = "all")
  
  expect_equal(
    round(raster::getValues(getSatDataLayer(t1, bcde = "B002n_RAD"))[50],4), 
    round(62.95932, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t2, bcde = "B002n_REF"))[50],4), 
    round( 0.1179185, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t3, bcde = "B002n_REF_NoSZEN"))[50],4), 
    round(0.10124, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t4, bcde = "B010n_BT"))[50],4), 
    round(303.5775, 4))
  expect_equal(
    round(raster::getValues(getSatDataLayer(t5, bcde = "B009n_RAD"))[50],4), 
    round(0.1295196, 4))
})

test_that("Depricated satCalib and calibLinear works as expected", {
  test_that("Depricated satCalib and calibLinear works as expected", {
    path <- system.file("extdata", 
                        package = "satellite")
    files <- list.files(path, 
                        pattern = glob2rx("LC8*.tif"), 
                        full.names = TRUE)
    sat <- satellite(files)
    bcde <- "B002n"
    
    t1 <- calibLinear(band = getSatDataLayer(sat, bcde),
                      mult = getSatREFM(sat, bcde),
                      add = getSatREFA(sat, bcde))
    
    t2 <- calibLinear(band = getSatDataLayer(sat, bcde),
                      mult = getSatREFM(sat, bcde),
                      add = getSatREFA(sat, bcde),
                      szen = getSatSZEN(sat, bcde))
    
    t3 <- calibLinear(band = getSatDataLayer(sat, bcde),
                      mult = getSatRADM(sat, bcde),
                      add = getSatRADA(sat, bcde))
    
    bcde <- "B010n"
    t4 <- calibLinear(band = getSatDataLayer(sat, bcde),
                      mult = getSatRADM(sat, bcde),
                      add = getSatRADA(sat, bcde),
                      k1 = getSatBTK1(sat, bcde),
                      k2 = getSatBTK2(sat, bcde))
    
    expect_equal(round(raster::getValues(t1)[50],3), round(0.10124,3))
    expect_equal(round(raster::getValues(t2)[50],3), round(0.1179185,3))
    expect_equal(round(raster::getValues(t3)[50],3), round(62.959,3))
    expect_equal(round(raster::getValues(t4)[50],3), round(303.5775,3))
  })
})