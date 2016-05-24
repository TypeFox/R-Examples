##' testthat tests for bad calls to function
##'
##' @author Raphael Winkelmann
context("test bad calls")

test_that("length(listOfFiles) > 1 and ToFile=F causes Error", {
  
  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)
  
  for (func in names(wrasspOutputInfos)){
    funcFormals = formals(func)
    funcFormals$listOfFiles = wavFiles
    funcFormals$ToFile = FALSE
    expect_that(do.call(func,as.list(funcFormals)), throws_error())
  }
  
})
