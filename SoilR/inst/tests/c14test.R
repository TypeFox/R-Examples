#
# vim:set ff=unix expandtab ts=2 sw=2:
#!/usr/bin/Rscript
source("prolog.R")
c14tests<- defineTestSuite(
   name="c14tests",
   #dirs = file.path(.path.package(package="SoilR"),
   #"tests"),
   dirs=".",
   testFileRegexp = "^runit.automatic.c14.+\\.[rR]",
   testFuncRegexp = "^test.+",
   rngKind = "Marsaglia-Multicarry",
   rngNormalKind = "Kinderman-Ramage"
)

testResult <- runTestSuite(c14tests)
printTextProtocol(testResult)
#produce exitstatus ne 0 for buildbot to notice
ef=getErrors(testResult)
n=ef$nErr+ef$nFail
if (n>0) {stop(1)}
