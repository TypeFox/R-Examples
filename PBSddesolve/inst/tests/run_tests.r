#Source this file to run all tests

require( "PBSddesolve" )
require( "RUnit" )

testsuite.dde <- defineTestSuite("dde",
	dirs = file.path(.path.package(package="PBSddesolve"), "tests"),
	testFileRegexp = "^runit.+\\.r",
	testFuncRegexp = "^test.+",
	rngKind = "Marsaglia-Multicarry",
	rngNormalKind = "Kinderman-Ramage")


testResult <- runTestSuite(testsuite.dde)
printTextProtocol(testResult)