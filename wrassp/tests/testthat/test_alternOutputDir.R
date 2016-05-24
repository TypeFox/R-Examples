##' testthat test for alternative output dir
##'
##' @author Raphael Winkelmann
context("test alternative output dir")

test_that("all signal processing functions run without errors on audio files", {
  
  altDir = tempdir()
  
  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)
  
  for (func in names(wrasspOutputInfos)){
    funcFormals = formals(func)
    funcFormals$listOfFiles = wavFiles
    funcFormals$outputDirectory = altDir
    res = do.call(func,as.list(funcFormals))
    expect_that(res, equals(NULL))
  }
  
})
