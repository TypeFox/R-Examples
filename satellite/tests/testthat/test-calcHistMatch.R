# devtools::test(".", "calcHistMatch")
context("calcHistMatch")


#-------------------------------------------------------------------------------
test_that("calcHistMatch for raster works as expected", {
  path <- system.file("extdata", package = "satellite")
  
  files <- list.files(path, pattern = glob2rx("LE7*.tif"), full.names = TRUE)
  l7 <- satellite(files)
  
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  l8 <- satellite(files)
  
  x <- getSatDataLayer(l7, "B002n")
  target <- getSatDataLayer(l8, "B003n")

#   Rsenal::histmatchRaster(x, target)
  
  
  })