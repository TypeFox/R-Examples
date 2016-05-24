#
# vim:set ff=unix expandtab ts=2 sw=2:
#!/usr/bin/Rscript
library("RUnit")
library("SoilR")
#source_tests <- defineTestSuite(
#   name="sourc_Tests",
#   dirs = file.path(.path.package(package="SoilR"), "tests"),
#   #dirs=".",
#   testFileRegexp = "^runit.+\\.[rR]",
#   testFuncRegexp = "^test.+",
#   rngKind = "Marsaglia-Multicarry",
#   rngNormalKind = "Kinderman-Ramage"
#)
#package_tests <- defineTestSuite(
#   name="package_Tests",
#   dirs = file.path(.path.package(package="SoilR"), "tests"),
#   #dirs=".",
#   testFileRegexp = "^package.+\\.[rR]",
#   testFuncRegexp = "^test.+",
#   rngKind = "Marsaglia-Multicarry",
#   rngNormalKind = "Kinderman-Ramage"
#)
#
#testResult <- runTestSuite(testSuites=list(source_tests,package_tests))
#printTextProtocol(testResult)
#
##produce exitstatus ne 0 for buildbot to notice
#ef=getErrors(testResult)
#n=ef$nErr+ef$nFail
#if (n>0) {stop(1)}
