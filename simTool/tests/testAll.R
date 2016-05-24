require(RUnit)
require(simTool)
testsuite.all <- defineTestSuite("ALL", dirs = ".", testFileRegexp = "^runit.+\\.R", testFuncRegexp = "^test.+")#,
#                                 rngKind = "Marsaglia-Multicarry",
#                                 rngNormalKind = "Kinderman-Ramage")
printTextProtocol(rts <- runTestSuite(testsuite.all))
printHTMLProtocol(rts, fileName="out.html")
