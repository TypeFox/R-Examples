# Run all unit test files
# 
# Author: modelf
###############################################################################

library("RUnit")
library(mcr)

options(warn=1)

testSuite <- defineTestSuite(name="MCReg",
		dirs=".",
		testFileRegexp="runit.*\\.R$",
		rngKind="default",
		rngNormalKind="default")

testData <- runTestSuite(testSuite, verbose=0L)
printTextProtocol(testData, showDetails=FALSE)
printHTMLProtocol(testData,file="testProtocol.html")
