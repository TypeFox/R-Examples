library(futureheatwaves)
context("Acquiring directory structure")

dataFolder <- system.file("extdata/cmip5", package = "futureheatwaves")
coordinateFilenames <- "latitude_longitude_NorthAmerica_12mo.csv"
tasFilenames <- "tas_NorthAmerica_12mo.csv"
timeFilenames <- "time_NorthAmerica_12mo.csv"

mods_1 <- "bcc1"
mods_2 <- "ccsm"
mods_3 <- c("bcc1", "ccsm")

finalList <- acquireDirectoryStructure(dataFolder = dataFolder,
                                       coordinateFilenames = coordinateFilenames,
                                       tasFilenames = tasFilenames,
                                       timeFilenames = timeFilenames,
                                       models_to_run = "all",
                                       dataDirectories = list(
                                               "historical" = c(1990, 1999),
                                               "rcp85" = c(2060, 2079)),
                                       threshold_ensemble = "r1i1p1",
                                       thresholdBoundaries = c(1990, 1999))

finalList_1 <- acquireDirectoryStructure(dataFolder = dataFolder,
                                       coordinateFilenames = coordinateFilenames,
                                       tasFilenames = tasFilenames,
                                       timeFilenames = timeFilenames,
                                       models_to_run = mods_1,
                                       dataDirectories = list(
                                               "historical" = c(1990, 1999),
                                               "rcp85" = c(2060, 2079)),
                                       threshold_ensemble = "r1i1p1",
                                       thresholdBoundaries = c(1990, 1999))

finalList_2 <- acquireDirectoryStructure(dataFolder = dataFolder,
                                         coordinateFilenames = coordinateFilenames,
                                         tasFilenames = tasFilenames,
                                         timeFilenames = timeFilenames,
                                         models_to_run = mods_2,
                                         dataDirectories = list(
                                                 "historical" = c(1990, 1999),
                                                 "rcp85" = c(2060, 2079)),
                                         threshold_ensemble = "r1i1p1",
                                         thresholdBoundaries = c(1990, 1999))

finalList_3 <- acquireDirectoryStructure(dataFolder = dataFolder,
                                         coordinateFilenames = coordinateFilenames,
                                         tasFilenames = tasFilenames,
                                         timeFilenames = timeFilenames,
                                         models_to_run = mods_3,
                                         dataDirectories = list(
                                                 "historical" = c(1990, 1999),
                                                 "rcp85" = c(2060, 2079)),
                                         threshold_ensemble = "r1i1p1",
                                         thresholdBoundaries = c(1990, 1999))

test_that("User can specify which models to run", {
        expect_equal(as.vector(sapply(finalList, FUN = function(x) x[[1]])),
                     c("bcc1", "ccsm"))
        expect_equal(as.vector(sapply(finalList_1, FUN = function(x) x[[1]])),
                     "bcc1")
})

finalList_4 <- acquireDirectoryStructure(dataFolder = dataFolder,
                                         coordinateFilenames = coordinateFilenames,
                                         tasFilenames = tasFilenames,
                                         timeFilenames = timeFilenames,
                                         models_to_run = "all",
                                         dataDirectories = list(
                                                 "historical" = c(1990, 1999),
                                                 "rcp85" = c(2060, 2079)),
                                         threshold_ensemble = "r1i1p1",
                                         thresholdBoundaries = c(1990, 1999))
finalList_5 <- acquireDirectoryStructure(dataFolder = dataFolder,
                                         coordinateFilenames = coordinateFilenames,
                                         tasFilenames = tasFilenames,
                                         timeFilenames = timeFilenames,
                                         models_to_run = "all",
                                         dataDirectories = list(
                                                 "historical" = c(1990, 1999),
                                                 "rcp85" = c(2060, 2079)),
                                         threshold_ensemble = "r2i1p1",
                                         thresholdBoundaries = c(1990, 1999))
finalList_6 <- acquireDirectoryStructure(dataFolder = dataFolder,
                                         coordinateFilenames = coordinateFilenames,
                                         tasFilenames = tasFilenames,
                                         timeFilenames = timeFilenames,
                                         models_to_run = "all",
                                         dataDirectories = list(
                                                 "historical" = c(1990, 1999),
                                                 "rcp85" = c(2060, 2079)),
                                         threshold_ensemble = "r2i1p1",
                                         thresholdBoundaries = c(2070, 2070))

test_that("acquireDirectoryStructure will weed out models without right ens. member", {
        expect_equal(as.vector(sapply(finalList_4, FUN = function(x) x[[1]])),
                     c("bcc1", "ccsm"))
        expect_equal(sapply(finalList_5, FUN = function(x) x[[1]]),
                     list())
        expect_equal(as.vector(sapply(finalList_6, FUN = function(x) x[[1]])),
                     "ccsm")
})
