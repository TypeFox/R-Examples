# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################

### Don't run if in CHECK mode
if(!grepl("SafeQuant\\.Rcheck",getwd())){
	
	setwd(dirname(sys.frame(1)$ofile))
	source("initTestSession.R")
	source("testParser.R")
	source("testIdentificationAnalysis.R")
	source("testExpressionAnalysis.R")


	source("testSafeQuantAnalysis.R")
	source("testTMT.R")

	source("testUserOptions.R")
	source("testGraphics.R")
	
	cat("\n ---------------------- ALL TESTS PASSED ------------------ \n")
	
}






