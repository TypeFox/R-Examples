#################################################
#
# This file provides a basic testing bed that 
# creates a window for each window description
# text file found in the same directory
#
#
# FOR AN INTERACTIVE DEMO USE:   testWidgets()
#
#
# To run: source("start_test.r")
#
#################################################

require("PBSmodelling")
x <- dir()
x <- x[grep("\\.txt$",x)]
for(i in 1:length(x)) {
	
	len<-nchar(x[i], "chars")+21
	for (j in 1:len) cat("-")
	cat("\n")
	cat("Creating window from "); cat(x[i]); cat("\n")
	for (j in 1:len) cat("-")
	cat("\n")
	
	cat("\n\n")
	
	tt<-createWin(x[i])
	print(getWinVal())

	cat("push <enter> to continue.\n")
	scan(nmax=1, what=character(), quiet=TRUE)

	closeWin()
}