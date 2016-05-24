
if(require("RUnit", quietly=TRUE)) {
  ## Run all the tests
  wd <- getwd()
  testsuite <- defineTestSuite("PCICt", dirs=wd, testFileRegexp = "^test_functions.R$", testFuncRegexp = "^PCICt.test.+")
  PCICt.test.result <- runTestSuite(testsuite, useOwnErrorHandler=F)
  printTextProtocol(PCICt.test.result)
}
