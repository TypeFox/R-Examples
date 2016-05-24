# devtools::test(".", "calcPathRadDOS")
context("calcPathRadDOS")


#-------------------------------------------------------------------------------
test_that("calcPathRadDOS for numeric works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)
  sat <- satTOAIrrad(sat, method = "Model")
  
  bcde <- "B002n"
  t1 <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
                       bnbr = getSatLNBR(sat, bcde),
                       band_wls = data.frame(LMIN = getSatLMIN(sat, getSatBCDESolar(sat)), 
                                             LMAX = getSatLMAX(sat, getSatBCDESolar(sat))),
                       radm = getSatRADM(sat, getSatBCDESolar(sat)),
                       rada = getSatRADA(sat, getSatBCDESolar(sat)),
                       szen = getSatSZEN(sat, getSatBCDESolar(sat)),
                       esun = getSatESUN(sat, getSatBCDESolar(sat)),
                       model = "DOS2",
                       scat_coef = -4)
  
  t2 <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
                       bnbr = getSatLNBR(sat, bcde),
                       band_wls = data.frame(LMIN = getSatLMIN(sat, getSatBCDESolar(sat)), 
                                             LMAX = getSatLMAX(sat, getSatBCDESolar(sat))),
                       radm = getSatRADM(sat, getSatBCDESolar(sat)),
                       rada = getSatRADA(sat, getSatBCDESolar(sat)),
                       szen = getSatSZEN(sat, getSatBCDESolar(sat)),
                       esun = getSatESUN(sat, getSatBCDESolar(sat)),
                       model = "DOS2",
                       scat_coef = -2)
  
  t3 <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
                       bnbr = getSatLNBR(sat, bcde),
                       band_wls = data.frame(LMIN = getSatLMIN(sat, getSatBCDESolar(sat)), 
                                             LMAX = getSatLMAX(sat, getSatBCDESolar(sat))),
                       radm = getSatRADM(sat, getSatBCDESolar(sat)),
                       rada = getSatRADA(sat, getSatBCDESolar(sat)),
                       szen = getSatSZEN(sat, getSatBCDESolar(sat)),
                       esun = getSatESUN(sat, getSatBCDESolar(sat)),
                       model = "DOS2",
                       scat_coef = -1)  
  
  t4 <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
                       bnbr = getSatLNBR(sat, bcde),
                       band_wls = data.frame(LMIN = getSatLMIN(sat, getSatBCDESolar(sat)), 
                                             LMAX = getSatLMAX(sat, getSatBCDESolar(sat))),
                       radm = getSatRADM(sat, getSatBCDESolar(sat)),
                       rada = getSatRADA(sat, getSatBCDESolar(sat)),
                       szen = getSatSZEN(sat, getSatBCDESolar(sat)),
                       esun = getSatESUN(sat, getSatBCDESolar(sat)),
                       model = "DOS2",
                       scat_coef = -0.7)
  
  t5 <- calcPathRadDOS(x = min(getValues(getSatDataLayer(sat, bcde))),
                       bnbr = getSatLNBR(sat, bcde),
                       band_wls = data.frame(LMIN = getSatLMIN(sat, getSatBCDESolar(sat)), 
                                             LMAX = getSatLMAX(sat, getSatBCDESolar(sat))),
                       radm = getSatRADM(sat, getSatBCDESolar(sat)),
                       rada = getSatRADA(sat, getSatBCDESolar(sat)),
                       szen = getSatSZEN(sat, getSatBCDESolar(sat)),
                       esun = getSatESUN(sat, getSatBCDESolar(sat)),
                       model = "DOS2",
                       scat_coef = -0.5)  
  
#   expect_equal(round(t1[1],3), c("B001n" = round(60.16885,3)))
#   expect_equal(round(t2[3],3), c("B003n" = round(29.51984,3)))
#   expect_equal(round(t3[4],3), c("B004n" = round(30.09144,3)))
#   expect_equal(round(t4[5],3), c("B005n" = round(28.29916,3)))
#   expect_equal(round(t5[6],3), c("B006n" = round(24.61562,3)))
#   
  #   c(coef-4, coef-2, coef-1  coef-0.7  coef-0.5
  #   1  59.8861832583 50.01006637 45.566175 44.293677 43.460494
  #   2  41.2138869229 41.21388692 41.213887 41.213887 41.213887
  #   3  20.1958822971 29.24265456 34.915089 36.792434 38.092498
  #   4   9.2349957403 20.77494563 29.856403 33.175119 35.564547
  #   5   1.9468918394 11.78153912 23.199239 28.156769 31.982744
  #   6  -0.2445099344  3.48253659 13.136381 19.160386 24.579342
  #   7  -0.1021297427  1.98553024  9.853833 15.682434 21.341225
  #   8  17.2188865978 26.78113766 33.418198 35.724149 37.352906
  #   9  -0.2151743576  4.74046853 15.244152 21.225972 26.397953
  #   10 -0.0002675506  0.08882900  2.030244  5.183168  9.680882
  #   11 -0.0001497598  0.07333347  1.843295  4.843989  9.223823
})


#-------------------------------------------------------------------------------
test_that("calcPathRadDOS for Satellite works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)
  
  sat_pathrad <- calcPathRadDOS(sat, model = "DOS2", esun_method = "RadRef")
  
#   expect_equal(round(getSatPRAD(sat_pathrad, bcde = "B002n"),3), 
#                round(c(B002n = 42.064), 3))
#   expect_equal(round(getSatPRAD(sat_pathrad, bcde = "B009n"),3), 
#                round(c(B009n = -0.185), 3))
})


#-------------------------------------------------------------------------------
  test_that("Deprecated satPathRadDOS for Satellite works as expected", {
    path <- system.file("extdata", package = "satellite")
    files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
    sat <- satellite(files)
    
    sat_pathrad <- satPathRadDOS(sat, atmos_model = "DOS2", 
                                 esun_mode = "RadRef")
    
#     expect_equal(round(getSatPRAD(sat_pathrad, bcde = "B002n"),3), 
#                  round(c(B002n = 42.064), 3))
#     expect_equal(round(getSatPRAD(sat_pathrad, bcde = "B009n"),3), 
#                  round(c(B009n = -0.185), 3))
  })