#
# vim:set ff=unix expandtab ts=2 sw=2:
library("RUnit")
c2f <- function(c) return(9/5 * c + 32)

test.c2f <- function() {
checkEquals(c2f(0), 32)
checkEquals(c2f(10), 50)
checkException(c2f("xx"))
}

testsuite.c2f <- defineTestSuite("c2f",
dirs = file.path(.path.package(package="RUnit"),
"examples"),
testFileRegexp = "^runit.+\\.r",
testFuncRegexp = "^test.+",
rngKind = "Marsaglia-Multicarry",
rngNormalKind = "Kinderman-Ramage")

testResult <- runTestSuite(testsuite.c2f)
printTextProtocol(testResult)

runTestFile(file.path(.path.package(package="RUnit"),
"examples/runitc2f.r"))

