#Source this file from the current directory to run all tests

if( require( "RUnit" ) ) {
	wd <- getwd()
	testsuite <- defineTestSuite("PBSmodelling",
		dirs = ".",
		testFileRegexp = "^runit.+\\.r$",
		testFuncRegexp = "^test.+",
		rngKind = "Marsaglia-Multicarry",
		rngNormalKind = "Kinderman-Ramage")
	
	testResult <- runTestSuite(testsuite)
	printTextProtocol(testResult)
	setwd( wd )
} else {
	print( "RUnit is not installed - cannot run unit tests" )
}