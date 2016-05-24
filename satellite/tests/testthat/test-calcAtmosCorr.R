# devtools::test(".", "calcAtmosCorr")
context("calcAtmosCorr")

#-------------------------------------------------------------------------------
test_that("calcAtmosCorr for RasterLayers works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)
  sat <- calcTOAIrradRadRef(sat)
  
  bcde <- "B002n"
  path_rad <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
                             bnbr = getSatLNBR(sat, bcde),
                             band_wls = 
                               data.frame(LMIN = 
                                            getSatLMIN(sat, 
                                                       getSatBCDESolar(sat)), 
                                          LMAX = 
                                            getSatLMAX(sat, 
                                                       getSatBCDESolar(sat))),
                             radm = getSatRADM(sat, getSatBCDESolar(sat)),
                             rada = getSatRADA(sat, getSatBCDESolar(sat)),
                             szen = getSatSZEN(sat, getSatBCDESolar(sat)),
                             esun = getSatESUN(sat, getSatBCDESolar(sat)),
                             model = "DOS2")
  
  sensor_rad <- convSC2Rad(x = getSatDataLayer(sat, bcde),
                                mult = getSatRADM(sat, bcde),
                                add = getSatRADA(sat, bcde))
  
  ref_atmos <- calcAtmosCorr(x = sensor_rad,
                             path_rad = path_rad[names(path_rad) == bcde],
                             esun = getSatESUN(sat, bcde),
                             szen = getSatSZEN(sat, bcde), 
                             model = "DOS2")
  
  expect_equal(round(raster::getValues(ref_atmos)[50],4), round(0.0455, 4))
})



#-------------------------------------------------------------------------------
test_that("calcAtmosCorr for Satellite works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)
  sat_atmos <- calcAtmosCorr(sat, model = "DOS2", esun_method = "RadRef")
  
  expect_equal(
    round(raster::getValues(
      getSatDataLayer(sat_atmos, bcde = "B002n_REF_AtmosCorr"))[50],4), 
    round(0.0441, 4))
})


#-------------------------------------------------------------------------------
test_that("Depricated satAtmosCorr for Satellite works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)
  sat_atmos <- satAtmosCorr(sat, atmos_model = "DOS2", esun_mode = "RadRef")
  
  expect_equal(
    round(raster::getValues(
      getSatDataLayer(sat_atmos, bcde = "B002n_REF_AtmosCorr"))[50],4), 
    round(0.0441, 4))
})