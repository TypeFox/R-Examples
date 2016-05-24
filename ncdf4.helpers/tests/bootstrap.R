if(require("RUnit", quietly=TRUE)) {
  library(ncdf4.helpers)
  library(proj4)

  ## Run all the tests
  wd <- getwd()
  testsuite <- defineTestSuite("ncdf4.helpers", dirs=wd, testFileRegexp = "^test_.+.R$", testFuncRegexp = "^ncdf4.helpers.test.+")
  ncdf4.helpers.test.result <- runTestSuite(testsuite, useOwnErrorHandler=F)
  printTextProtocol(ncdf4.helpers.test.result)
  stopifnot(ncdf4.helpers.test.result$ncdf4.helpers$nFail == 0 && ncdf4.helpers.test.result$ncdf4.helpers$nErr == 0)    
}
