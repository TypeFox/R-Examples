## when the tests codes have not published yet
library(RUnit)
library(VarianceGamma)
data(vgParam)
testParam <- vgSmallShape
testParam <- vgTinyNu
#####--------------------------------------------------------------------####
## run level 2 tests for all functions
TestSuitevgL2 <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = "E:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel2Tests.R")
testResultL2 <- runTestSuite(TestSuitevgL2)
## prints detailed text protocol
printTextProtocol(testResultL2, showDetails = TRUE)


#####--------------------------------------------------------------------####

## run level 2 tests for each functions group
### dpqr group - dvg, ddvg, pvg, qvg, rvg, vgCalcRange, vgBreaks
TestSuitevgL2dpqr <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = "E:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel2Tests.R", testFuncRegexp = "test.vgL2dpqr")
testResultL2dpqr <- runTestSuite(TestSuitevgL2dpqr)
## prints detailed text protocol
printTextProtocol(testResultL2dpqr, showDetails = TRUE)


### moments group - vgMode, vgMean, vgVar, vgSkew, vgKurt, vgMom
TestSuitevgL2moments <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = "C:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel2Tests.R", testFuncRegexp = "test.vgL2moments")
testResultL2moments <- runTestSuite(TestSuitevgL2moments)
## prints detailed text protocol
printTextProtocol(testResultL2moments, showDetails = TRUE)


### fitting group - vgFitStartMoM, vgFitStart, vgFit
TestSuitevgL2fitting <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = "C:/vg/vg0.2-2/unitTests/testFunction",
  testFileRegexp = "runit.vgLevel2Tests.R", testFuncRegexp = "test.vgL2fitting")
testResultL2fitting <- runTestSuite(TestSuitevgL2fitting)
## prints detailed text protocol
printTextProtocol(testResultL2fitting, showDetails = TRUE)



#####--------------------------------------------------------------------####

##when the test files and the R codes already published
library(RUnit)
library(VarianceGamma)
data(vgParam)
testParam <- vgLargeParam
vgTestSuiteL2 <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
testFileRegexp = "runit.vglevel2Tests.R")
testResult <- runTestSuite(vgTestSuiteL2)
## prints detailed text protocol
printTextProtocol(testResult, showDetails = TRUE)

## run level 2 tests for each functions group
### dpqr group - dvg, ddvg, pvg, qvg, rvg, vgCalcRange, vgBreaks
TestSuitevgL2dpqr <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel2Tests.R", testFuncRegexp = "test.vgL2dpqr")
testResultL2dpqr <- runTestSuite(TestSuitevgL2dpqr)
## prints detailed text protocol
printTextProtocol(testResultL2dpqr, showDetails = TRUE)


### moments group - vgMode, vgMean, vgVar, vgSkew, vgKurt, vgMom
TestSuitevgL2moments <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel2Tests.R", testFuncRegexp = "test.vgL2moments")
testResultL2moments <- runTestSuite(TestSuitevgL2moments)
## prints detailed text protocol
printTextProtocol(testResultL2moments, showDetails = TRUE)


### fitting group - vgFitStartMoM, vgFitStart, vgFit
TestSuitevgL2fitting <- defineTestSuite(name = "VG Level 2 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel2Tests.R", testFuncRegexp = "test.vgL2fitting")
testResultL2fitting <- runTestSuite(TestSuitevgL2fitting)
## prints detailed text protocol
printTextProtocol(testResultL2fitting, showDetails = TRUE)