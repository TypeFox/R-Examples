
runBLCOPTests <- function(testPath = BLCOPOptions("unitTestPath"), protocolFile = "BLCOPTests.html", writeProtocol = FALSE)
{
	BLTestSuite <- defineTestSuite(name = "Black-Litterman / COP unit tests", dirs = testPath) 
    testResults <- runTestSuite(BLTestSuite)
    if(writeProtocol)
		printHTMLProtocol(testResults, fileName = protocolFile)
	testResults
}