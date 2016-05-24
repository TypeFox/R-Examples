# Run all unit tests, i.e. all checks of all test-functions
# 
# Author: Andre Schützenmeister
###############################################################################

library(VCA)
library(RUnit)

# enable/disable time-consuming model 1 and model 2 testcases with real-world data
# see file "runit.VCAinference.R"

realWorldModel2 <- FALSE
realWorldModel1 <- FALSE

options(warn=1)

# test function regexpr fits to string "TFxyz" which are used as identifiers for easier referencing

testSuite <- defineTestSuite(name="VCA", dirs=".",
                             testFileRegexp="runit.*\\.R$",
							 testFuncRegexp = "^TF[[:digit:]]{3}.+",					# use custom regexpr for test functions
                             rngKind="default",
                             rngNormalKind="default")
                     
testData <- runTestSuite(testSuite, verbose=0L)

sInfo <- sessionInfo()
cat("Test Summary Report R-Paket VCA", paste("V", sInfo$otherPkgs[["VCA"]]$Version, sep=""), file="./VCA_UnitTest_Protocol.txt", append=FALSE)
cat("\n-------------------------------------", file="./VCA_UnitTest_Protocol.txt", append=TRUE )
cat("\n\n\n1) Package Description:", file="./VCA_UnitTest_Protocol.txt", append=TRUE)
cat("\n-----------------------\n\n", file="./VCA_UnitTest_Protocol.txt", append=TRUE)
capture.output(print(sInfo$otherPkgs[["VCA"]]), file="./VCA_UnitTest_Protocol.txt", append=TRUE)
cat("\n\n\n2) Test Environment:", file="./VCA_UnitTest_Protocol.txt", append=TRUE)
cat("\n--------------------\n", file="./VCA_UnitTest_Protocol.txt", append=TRUE)
sinfo <- Sys.info()
snam <- names(sinfo)
for(i in 1:length(sinfo))
{
	cat(paste("\n", snam[i], paste(rep(" ", 20-nchar(snam[i])), collapse=""), sep=""),":\t", sinfo[i], file="./VCA_UnitTest_Protocol.txt", append=TRUE)
}

cat("\n\n\n\n3) Test Protocol:", file="./VCA_UnitTest_Protocol.txt", append=TRUE)
cat("\n-----------------\n\n", file="./VCA_UnitTest_Protocol.txt", append=TRUE)

printTextProtocol(testData, showDetails=FALSE)
capture.output(printTextProtocol(testData, showDetails=TRUE), file="./VCA_UnitTest_Protocol.txt", append=TRUE)
printHTMLProtocol(testData, file="./VCA_UnitTest_Protocol.html")

options(warn=0)

