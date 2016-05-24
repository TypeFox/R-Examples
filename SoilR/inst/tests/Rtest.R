#!/usr/bin/Rscript
# vim:set ff=unix expandtab ts=2 sw=2:
source("prolog.R")
tfr  <- "^runit\\..*\\.R"
#fl <- list.files(pattern=tfr)
#for (fn in fl){
#  print(fn)
#  source(fn)
#}
alltests <- defineTestSuite(
   name="allTests",
   dirs=c(".","protected"),
   testFileRegexp = tfr,
   
   #testFuncRegexp = "^test.TwopSerial_linear_vs_nonlinear"
   #"^test.FourpSerial_1"
   #"test.TwopParallel_ZeroInput"
   #"^test.TwopFeedback"
   #"^test.TimeMapInterface"
   #"^test.LowVerticalRatesPaper" 
   #"^test.check.pass"
   #"test.ModelInit"
   #"ptest.ModelOperators"
   #"test.ParallelModel"
   #"test.TwopSerial_linear_vs_nonlinear"
   #"test.SoilRPaper1"
   #"test.FourpSerial_1"
   #"test.BoundFc"
   #"test.ThreepFeedbackModel14|test.ThreepParallelModel14|test.ThreepSeriesModel14|test.TwopFeedbackModel14|test.TwopParallelModel14|test.TwopSeriesModel14"
   #"test.LowVerticalRatesPaper|test.ModelInit|test.SoilRPaper1"
   "test.LowVerticalRatesPaper|test.SoilRPaper1"
   #"test.Deprecation"
   #"test.GaudinskiModel14"
   #"test.MC"
)

testResult <- runTestSuite(alltests)
printTextProtocol(testResult)
#produce exitstatus ne 0 for buildbot to notice
ef=getErrors(testResult)
n=ef$nErr+ef$nFail
if (n>0) {stop(1)}
