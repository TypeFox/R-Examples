### R code from vignette source 'HistogramTools-UnitTests.Rnw'

###################################################
### code chunk number 1: HistogramTools-UnitTests.Rnw:25-30
###################################################
options(width=50)
library(HistogramTools)
library(RUnit)
set.seed(0)
ht.version <- packageDescription("HistogramTools")$Version


###################################################
### code chunk number 2: HistogramTools-UnitTests.Rnw:41-49
###################################################
# Define tests
testSuite <- defineTestSuite(
  name="HistogramTools Unit Tests",
  dirs=system.file("unitTests", package = "HistogramTools"),
  testFuncRegexp = "^[Tt]est.+")

# Run tests
tests <- runTestSuite(testSuite)


###################################################
### code chunk number 3: HistogramTools-UnitTests.Rnw:52-54
###################################################
# Print results
printTextProtocol(tests)


