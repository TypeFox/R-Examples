
library(RUnit)
library(RSAP)
 
test.suite <- defineTestSuite("RSAP.core.tests",
                              dirs = file.path("tests"),
                              testFileRegexp = '^\\d+\\.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)

