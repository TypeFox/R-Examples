##
## runner.r - run RUnit tests included in package
##
## Based on code taken from the fBasics package.
##

pkg <- "fftw"

if (require("RUnit", quietly = TRUE)) {
  library(package=pkg, character.only = TRUE)
  if (!(exists("path") && file.exists(path)))
    path <- system.file("unitTests", package = pkg)
  
  ## --- Testing ---  
  ## Define tests
  testSuite <- defineTestSuite(name = paste(pkg, "unit testing"),
                               dirs = path)
  
  if (interactive()) {
    cat("Now have RUnit Test Suite 'testSuite' for package '", pkg, "' :\n", sep='')
    str(testSuite)
    cat('', "Consider doing",
        "\t  tests <- runTestSuite(testSuite)", "\nand later",
        "\t  printTextProtocol(tests)", '', sep = "\n")
  } else {
    tests <- runTestSuite(testSuite)
    
    if (file.access(path, 02) != 0) {
      tdir <- tempfile(paste(pkg, "unitTests", sep="_"))
      dir.create(tdir)
      pathReport <- file.path(tdir, "report")
      cat("RUnit reports saved to ", tdir, "/report.(txt|html)", sep = "")
    } else {
      pathReport <- file.path(path, "report")
    }
    
    ## Print Results:
    printTextProtocol(tests, showDetails = FALSE)
    printTextProtocol(tests, showDetails = FALSE,
                      fileName = paste(pathReport, "Summary.txt", sep = ""))
    printTextProtocol(tests, showDetails = TRUE,
                      fileName = paste(pathReport, ".txt", sep = ""))
    
    ## stop() if there are any failures i.e. FALSE to unit test.
    ## This will cause R CMD check to return error and stop
    tmp <- getErrors(tests)
    if (tmp$nFail > 0 | tmp$nErr > 0) {
      stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
                 ", R errors: ",  tmp$nErr, ")\n\n", sep=""))
    }
  }
} else {
  message("Unable to load R package 'RUnit'. No unit tests run.")
}
