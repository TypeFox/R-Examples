# ===========================================================================
# File: "runAllTests.R"
#                        Created: 2012-11-01 18:17:51
#              Last modification: 2012-11-07 14:46:36
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================


pkg <- "clusterCrit"

if (require("RUnit", quietly = TRUE)) {
    library(package=pkg, character.only = TRUE)
	
	unitTestDir <- system.file("unitTests", package = pkg)
 	clusTestDirs <- file.path(unitTestDir,c("intCriteria","extCriteria","others"))

    # Define tests
    testSuite <- defineTestSuite(name = paste(pkg, "unit testing"),
                                 dirs = clusTestDirs, 
								 testFuncRegexp = "^test\\..+")

	# Run the test suite
	tests <- runTestSuite(testSuite)

	# Set the output directory
	if (file.access(unitTestDir, 02) != 0) {
		# No write permission, create temporary directory
		reportsDir <- tempfile(paste(pkg, "unitTests", sep="_"))
		dir.create(reportsDir)
	} else {
		reportsDir <- unitTestDir
	}

	# Print the results
	printTextProtocol(tests, showDetails = FALSE)
	# Text version
	printTextProtocol(tests, showDetails = TRUE,
					  fileName = file.path(reportsDir, "report.txt"))
	# HTML version
	printHTMLProtocol(tests, fileName = file.path(reportsDir, "report.html"))

	# stop() if there are any failures i.e. FALSE to unit test.
	# This will cause R CMD check to return error and stop
	tmp <- getErrors(tests)
	if (tmp$nFail > 0 | tmp$nErr > 0) {
		stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
				   ", R errors: ",  tmp$nErr, ")\n\n", sep=""))
	} else {
		cat("RUnit reports are written in '", reportsDir, "':\n", sep = "")
		cat("\treport.txt\n")
		cat("\treport.html\n")
	}
} else {
    cat("R package 'RUnit' cannot be loaded -- no unit tests run\n",
        "for package", pkg,"\n")
}



