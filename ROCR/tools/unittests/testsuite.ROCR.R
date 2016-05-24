library(RUnit)
library(ROCR)

myTestSuite <- defineTestSuite("ROCR test suite", "unittests")
isValidTestSuite(myTestSuite)
testData <- runTestSuite(myTestSuite)

printTextProtocol(testData, showDetails=TRUE)
printHTMLProtocol(testData, "unittests/testresults.html")
