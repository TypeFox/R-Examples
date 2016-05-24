
## when the tests and the new R codes have not published yet
library(RUnit)
library(VarianceGamma)

#####--------------------------------------------------------------------####


## run level 1 tests for all functions
TestSuitevgL1 <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = "C:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel1Tests.R")
testResultL1 <- runTestSuite(TestSuitevgL1)
## prints detailed text protocol
printTextProtocol(testResultL1, showDetails = TRUE)


#####--------------------------------------------------------------------####
## run level 1 tests for each functions group

### dpqr group - dvg, ddvg, pvg, qvg, rvg, vgCalcRange, vgBreaks
TestSuitevgL1dpqr <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1dpqr")
testResultL1dpqr <- runTestSuite(TestSuitevgL1dpqr)
## prints detailed text protocol
printTextProtocol(testResultL1dpqr, showDetails = TRUE)


### moments group - vgMode, vgMean, vgVar, vgSkew, vgKurt, vgMom
TestSuitevgL1moments <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1moments")
testResultL1moments <- runTestSuite(TestSuitevgL1moments)
## prints detailed text protocol
printTextProtocol(testResultL1moments, showDetails = TRUE)


### fitting group - vgFitStartMoM, vgFitStart, vgFit
TestSuitevgL1fitting <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1fitting")
testResultL1fitting <- runTestSuite(TestSuitevgL1fitting)
## prints detailed text protocol
printTextProtocol(testResultL1fitting, showDetails = TRUE)

### parameter group - vgChangePars
TestSuitevgL1parameter <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1parameter")
testResultL1parameter <- runTestSuite(TestSuitevgL1parameter)
## prints detailed text protocol
printTextProtocol(testResultL1parameter, showDetails = TRUE)

#####--------------------------------------------------------------------####



##when the test files and the R codes already published
library(RUnit)
library(VarianceGamma)
vgTestSuiteL1 <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vglevel1Tests.R")
testResult <- runTestSuite(vgTestSuiteL1)
## prints detailed text protocol
printTextProtocol(testResult, showDetails = TRUE)

#####--------------------------------------------------------------------####
## run level 1 tests for each functions group

### dpqr group - dvg, ddvg, pvg, qvg, rvg, vgCalcRange, vgBreaks
TestSuitevgL1dpqr <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1dpqr")
testResultL1dpqr <- runTestSuite(TestSuitevgL1dpqr)
## prints detailed text protocol
printTextProtocol(testResultL1dpqr, showDetails = TRUE)


### moments group - vgMode, vgMean, vgVar, vgSkew, vgKurt, vgMom
TestSuitevgL1moments <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1moments")
testResultL1moments <- runTestSuite(TestSuitevgL1moments)
## prints detailed text protocol
printTextProtocol(testResultL1moments, showDetails = TRUE)


### fitting group - vgFitStartMoM, vgFitStart, vgFit
TestSuitevgL1fitting <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1fitting")
testResultL1fitting <- runTestSuite(TestSuitevgL1fitting)
## prints detailed text protocol
printTextProtocol(testResultL1fitting, showDetails = TRUE)

### parameter group - vgChangePars
TestSuitevgL1parameter <- defineTestSuite(name = "VG Level 1 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel1Tests.R", testFuncRegexp = "test.vgL1parameter")
testResultL1parameter <- runTestSuite(TestSuitevgL1parameter)
## prints detailed text protocol
printTextProtocol(testResultL1parameter, showDetails = TRUE)