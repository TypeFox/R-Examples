context("maskInvarFeatures")

test_that("maskInvarFeatures works as expected", {
  path <- system.file("extdata", 
                      package = "satellite")
  files <- list.files(path, 
                      pattern = glob2rx("LC8*.tif"), 
                      full.names = TRUE)
  sat <- satellite(files)

  t <- maskInvarFeatures(x = getSatDataLayer(sat, "B004n"), 
                    nir = getSatDataLayer(sat, "B005n"), 
                    swir = getSatDataLayer(sat, "B007n"))
})






