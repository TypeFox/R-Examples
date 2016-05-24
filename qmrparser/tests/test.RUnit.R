## unit tests will not be done if RUnit is not available

if(require("RUnit", quietly=TRUE)) {
  ## 
  ## --- Setup ---
  ##
  pkg <- "qmrparser" # <-- Change to package name!
                                        #3
  ## When in developing
  path <-  file.path(getwd(), "runit")
  if ( file.exists(path) ) {
    ## When in developing, developing environment is prioritised
    .libPaths(c(file.path(getwd(), ".."),.libPaths()))
  } else {
    ## Otherwise, check if package already installed
    path <- system.file(package=pkg, "tests", "runit")    
  }    
  ##
  ##
  ##
  cat("\nRunning unit tests\n")
  print(list(pkg=pkg, getwd=getwd(), pathToUnitTests=path))
  ##
  require(package=pkg, character.only=TRUE)
  ##
  wd <- setwd(path) ; path <- "." 
 
  ## If desired, load the name space to allow testing of private functions
  ## if (is.element(pkg, loadedNamespaces()))
  ##     attach(loadNamespace(pkg), name=paste("namespace", pkg, sep=":"), pos=3)
  ##
  ## or simply call PKG:::myPrivateFunction() in tests
 
  ## --- Testing ---
 
  ## Define tests
  testSuite <- defineTestSuite(
                name=paste(pkg, "unit testing"),
                dirs=path, 
                testFileRegexp = "^runit.+\\.[rR]$",
                testFuncRegexp = "^test.+")
  ## Run
  tests <- runTestSuite(testSuite)
 
  ## Default report name
  pathReport <- file.path(path, "report")
  dir.create(path=pathReport, recursive = TRUE)
 
  ## Report to stdout and text files
  cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
  printTextProtocol(tests, showDetails=TRUE,
                    fileName=file.path(pathReport, "report.txt"))
  printTextProtocol(tests, showDetails=FALSE,
                    fileName=file.path(pathReport, "Summary.txt"))
 
  ## Report to HTML file
  testFileToSFLinkMap <- function(testFileName, testDir = "..") {
    return(file.path("..",basename(testFileName)))
  }

  printHTMLProtocol(tests, fileName=file.path(pathReport, "report.html"),
                    testFileToLinkMap = testFileToSFLinkMap)
  
  setwd(wd)
  ## Return stop() to cause R CMD check stop in case of
  ##  - failures i.e. FALSE to unit tests or
  ##  - errors i.e. R errors
  tmp <- getErrors(tests)
  if(tmp$nFail > 0 | tmp$nErr > 0) {
    stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
               ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
  }
} else {
  warning("cannot run unit tests -- package RUnit is not available")
}
