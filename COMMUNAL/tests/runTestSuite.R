if(require(RUnit, quietly = TRUE)) {
  library(COMMUNAL)
  #runTestFile(file.path(.path.package(package="COMMUNAL"), "inst/tests/unittest.R"))
  testSuite <- defineTestSuite("COMMUNAL-test", dirs=file.path(path.package(package="COMMUNAL"), "tests"),
                               testFileRegexp = "^unittest.R", testFuncRegexp = "^test.+")
  testResult <- runTestSuite(testSuite)
  printTextProtocol(testResult)
}
