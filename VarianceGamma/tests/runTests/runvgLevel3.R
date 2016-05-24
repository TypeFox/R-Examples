## when the tests codes have not published yet
library(RUnit)
library(VarianceGamma)
data(vgParam)
testParam <- vgSmallShape
n <- 10000
nqp <- 100
N <- 100
Nfit <- 100
errorThresholdM <- 0.1
errorThresholdV <- 0.001
errorThresholdS <- 0.001
errorThresholdK <- 0.001
errorThresholdMom <- 0.01
errorThresholddpqrI <- 0.001
thresholddpqrR <- 0.001
errorThresholdFit <- 0.5
#####--------------------------------------------------------------------####
## run level 3 tests for all functions
TestSuitevgL3 <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/tests",
  testFileRegexp = "runit.vgLevel3Tests.R")
testResultL3 <- runTestSuite(TestSuitevgL3)
## prints detailed text protocol
printTextProtocol(testResultL3, showDetails = TRUE)


#####--------------------------------------------------------------------####

## run level 3 tests for each functions group
### dpqr group - dvg, ddvg, pvg, qvg, rvg, vgCalcRange, vgBreaks
TestSuitevgL3dpqrInversionqp <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3dpqr.R", testFuncRegexp = "test.vgL3dpqrInversionqp")
testResultL3dpqrInversionqp <- runTestSuite(TestSuitevgL3dpqrInversionqp)
## prints detailed text protocol
printTextProtocol(testResultL3dpqrInversionqp, showDetails = TRUE)

TestSuitevgL3dpqrInversionpq <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3dpqr.R", testFuncRegexp = "test.vgL3dpqrInversionpq")
testResultL3dpqrInversionpq <- runTestSuite(TestSuitevgL3dpqrInversionpq)
## prints detailed text protocol
printTextProtocol(testResultL3dpqrInversionpq, showDetails = TRUE)

TestSuitevgL3dpqrRandom <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3dpqr.R", testFuncRegexp = "test.vgL3dpqrRandom")
testResultL3dpqrRandom <- runTestSuite(TestSuitevgL3dpqrRandom)
## prints detailed text protocol
printTextProtocol(testResultL3dpqrRandom, showDetails = TRUE)


### moments group - vgMode, vgMean, vgVar, vgSkew, vgKurt, vgMom
TestSuitevgL3mean <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Mean")
testResultL3mean <- runTestSuite(TestSuitevgL3mean)
## prints detailed text protocol
printTextProtocol(testResultL3mean, showDetails = TRUE)

TestSuitevgL3var <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "F:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Var")
testResultL3var <- runTestSuite(TestSuitevgL3var)
## prints detailed text protocol
printTextProtocol(testResultL3var, showDetails = TRUE)

TestSuitevgL3skew <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "F:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Skew")
testResultL3skew <- runTestSuite(TestSuitevgL3skew)
## prints detailed text protocol
printTextProtocol(testResultL3skew, showDetails = TRUE)

TestSuitevgL3kurt <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "F:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Kurt")
testResultL3kurt <- runTestSuite(TestSuitevgL3kurt)
## prints detailed text protocol
printTextProtocol(testResultL3kurt, showDetails = TRUE)

TestSuitevgL3mom <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Mom")
testResultL3mom <- runTestSuite(TestSuitevgL3mom)
## prints detailed text protocol
printTextProtocol(testResultL3mom, showDetails = TRUE)


### fitting group - vgFitStartMoM, vgFitStart, vgFit
TestSuitevgL3fittingNM <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3fit.R", testFuncRegexp = "test.vgL3fitNM")
testResultL3fittingNM <- runTestSuite(TestSuitevgL3fittingNM)
## prints detailed text protocol
printTextProtocol(testResultL3fittingNM, showDetails = TRUE)

TestSuitevgL3fittingnlm <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3fit.R", testFuncRegexp = "test.vgL3fitnlm")
testResultL3fittingnlm <- runTestSuite(TestSuitevgL3fittingnlm)
## prints detailed text protocol
printTextProtocol(testResultL3fittingnlm, showDetails = TRUE)

TestSuitevgL3fittingBFGS <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/testFunction/runit.vgLevel3Tests",
  testFileRegexp = "runit.vgLevel3fit.R", testFuncRegexp = "test.vgL3fitBFGS")
testResultL3fittingBFGS <- runTestSuite(TestSuitevgL3fittingBFGS)
## prints detailed text protocol
printTextProtocol(testResultL3fittingBFGS, showDetails = TRUE)


#####--------------------------------------------------------------------####


##when the test files and the R codes already published
library(RUnit)
library(VarianceGamma)
data(vgParam)
testParam <- vgSmallShape
n <- 10000
nqp <- 100
N <- 100
Nfit <- 100
errorThresholdM <- 0.1
errorThresholdV <- 0.001
errorThresholdS <- 0.001
errorThresholdK <- 0.001
errorThresholdMom <- 0.01
errorThresholddpqrI <- 0.001
thresholddpqrR <- 0.001
errorThresholdFit <- 0.5
## run level 3 tests for all functions
TestSuitevgL3 <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = "G:/vg/vg0.2-2/unitTests/tests",
  testFileRegexp = "runit.vgLevel3Tests.R")
testResultL3 <- runTestSuite(TestSuitevgL3)
## prints detailed text protocol
printTextProtocol(testResultL3, showDetails = TRUE)


#####--------------------------------------------------------------------####

## run level 3 tests for each functions group
### dpqr group - dvg, ddvg, pvg, qvg, rvg, vgCalcRange, vgBreaks
TestSuitevgL3dpqrInversionqp <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3dpqr.R", testFuncRegexp = "test.vgL3dpqrInversionqp")
testResultL3dpqrInversionqp <- runTestSuite(TestSuitevgL3dpqrInversionqp)
## prints detailed text protocol
printTextProtocol(testResultL3dpqrInversionqp, showDetails = TRUE)

TestSuitevgL3dpqrInversionpq <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3dpqr.R", testFuncRegexp = "test.vgL3dpqrInversionpq")
testResultL3dpqrInversionpq <- runTestSuite(TestSuitevgL3dpqrInversionpq)
## prints detailed text protocol
printTextProtocol(testResultL3dpqrInversionpq, showDetails = TRUE)

TestSuitevgL3dpqrRandom <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3dpqr.R", testFuncRegexp = "test.vgL3dpqrRandom")
testResultL3dpqrRandom <- runTestSuite(TestSuitevgL3dpqrRandom)
## prints detailed text protocol
printTextProtocol(testResultL3dpqrRandom, showDetails = TRUE)


### moments group - vgMode, vgMean, vgVar, vgSkew, vgKurt, vgMom
TestSuitevgL3mean <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Mean")
testResultL3mean <- runTestSuite(TestSuitevgL3mean)
## prints detailed text protocol
printTextProtocol(testResultL3mean, showDetails = TRUE)

TestSuitevgL3var <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Var")
testResultL3var <- runTestSuite(TestSuitevgL3var)
## prints detailed text protocol
printTextProtocol(testResultL3var, showDetails = TRUE)

TestSuitevgL3skew <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Skew")
testResultL3skew <- runTestSuite(TestSuitevgL3skew)
## prints detailed text protocol
printTextProtocol(testResultL3skew, showDetails = TRUE)

TestSuitevgL3kurt <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Kurt")
testResultL3kurt <- runTestSuite(TestSuitevgL3kurt)
## prints detailed text protocol
printTextProtocol(testResultL3kurt, showDetails = TRUE)

TestSuitevgL3mom <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3moments.R", testFuncRegexp = "test.vgL3Mom")
testResultL3mom <- runTestSuite(TestSuitevgL3mom)
## prints detailed text protocol
printTextProtocol(testResultL3mom, showDetails = TRUE)


### fitting group - vgFitStartMoM, vgFitStart, vgFit
TestSuitevgL3fittingNM <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3fit.R", testFuncRegexp = "test.vgL3fitNM")
testResultL3fittingNM <- runTestSuite(TestSuitevgL3fittingNM)
## prints detailed text protocol
printTextProtocol(testResultL3fittingNM, showDetails = TRUE)

TestSuitevgL3fittingnlm <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3fit.R", testFuncRegexp = "test.vgL3fitnlm")
testResultL3fittingnlm <- runTestSuite(TestSuitevgL3fittingnlm)
## prints detailed text protocol
printTextProtocol(testResultL3fittingnlm, showDetails = TRUE)

TestSuitevgL3fittingBFGS <- defineTestSuite(name = "VG Level 3 RUnit Tests",
  dirs = system.file("testFunction", package = "VarianceGamma"),
  testFileRegexp = "runit.vgLevel3fit.R", testFuncRegexp = "test.vgL3fitBFGS")
testResultL3fittingBFGS <- runTestSuite(TestSuitevgL3fittingBFGS)
## prints detailed text protocol
printTextProtocol(testResultL3fittingBFGS, showDetails = TRUE)
