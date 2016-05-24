library("RUnit")
library("ff")
library("cluster")
library("golubEsets")
library("boot")
library("ShortRead")
library("e1071")
library("multtest")
library("stringdist")

# Automatically runs all R files under the unitTests directory. 
# Must make sure all necessary libraries are added to the top of this file when writing new tests.

filename <- paste("all_results",Sys.Date(),Sys.info()["sysname"],".log",sep = "_")
logFile <- file(filename)
sink(file = logFile, append = TRUE, type = c("output"), split = FALSE)
sink(file = logFile, append = TRUE, type = c("message"), split = FALSE)

allDirs  <- vector()
for (d in list.files("../unitTests")){
	thisDir = paste("../unitTests/",d, sep="")
	allDirs <- append(allDirs, thisDir, after = length(allDirs))
	for (nm in list.files(thisDir, pattern = "\\.[Rr]$")){	source(file.path(thisDir, nm))
	}
}
print(allDirs)
library("sprint")

print("*** System info ***")
sessionInfo()
R.version
Sys.info()
system("mpiexec -version")
print("*** End of system info ***")

# test.suite <- defineTestSuite("allUnitTests", dirs = file.path("../unitTests/ppam/"),testFileRegexp = '*.R')
test.suite <- defineTestSuite("allUnitTests", dirs = allDirs, testFileRegexp = '*.R')


# === Set up finished ===

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)

printHTMLProtocol(test.result, fileName = "All-unit-test-log.html")

print("*** End of tests ***")

pterminate()
quit()


