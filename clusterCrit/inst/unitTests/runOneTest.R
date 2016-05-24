# ===========================================================================
# File: "runAllTests.R"
#                        Created: 2012-11-01 18:17:51
#              Last modification: 2013-04-20 08:38:17
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================

testName <- commandArgs(trailingOnly = TRUE)
if (testName == "") {
	stop("Empty test name\n\n")
}

pkg <- "clusterCrit"

if (require("RUnit", quietly = TRUE)) {
    library(package=pkg, character.only = TRUE)
	
	# Find the criterion name
	testName <- paste("runit_",testName,".R",sep="")
	subdirs <- c("intCriteria","extCriteria","others")
	unitTestDir <- system.file("unitTests", package=pkg)
	found <- FALSE
	for (sub in subdirs) {
		testDir <- file.path(unitTestDir,sub)
		testList <- dir(path=testDir)
		if (testName %in% testList) {
			found <- TRUE
			break
		}
	}
	if (!found) {
		stop("\n\t\t>> can't find test for ",testName)
	}
		
	# Run the test file
	tst <- runTestFile(file.path(testDir,testName), 
		testFuncRegexp = "^test\\..+", 
		verbose = TRUE)

	# Set the output directory
	if (file.access(testDir, 02) != 0) {
		# No write permission, create temporary directory
		tdir <- tempfile(paste(pkg, "unitTests", sep="_"))
		reportsDir <- file.path(tdir, "reports")
	} else {
		unitTestDir <- system.file("unitTests", package=pkg)
		reportsDir <- file.path(unitTestDir, "reports")
	}
	dir.create(reportsDir)
	cat("RUnit reports are written into", reportsDir,"\n")

	# Print the results
	printTextProtocol(tst, showDetails = FALSE)
	# Text version
	printTextProtocol(tst, showDetails = TRUE,
					  fileName = file.path(reportsDir,paste("report_", testName, ".txt", sep = "")))
	# HTML version
	printHTMLProtocol(tst, fileName = file.path(reportsDir,paste("report_", testName, ".txt", sep = "")))

	# Stop in case of failures
	tmp <- getErrors(tst)
	if (tmp$nFail > 0 | tmp$nErr > 0) {
		stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
				   ", R errors: ",  tmp$nErr, ")\n\n", sep=""))
	}
    
} else {
    cat("R package 'RUnit' cannot be loaded -- no unit tests run\n",
        "for package", pkg,"\n")
}


