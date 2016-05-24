# ===========================================================================
# File: "runAllTests.R"
#                        Created: 2012-11-01 18:17:51
#              Last modification: 2012-11-02 13:41:59
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================

testName <- commandArgs(trailingOnly = TRUE)
if (testName == "") {
	stop("Empty test name\n\n")
}


## 
 # ------------------------------------------------------------------------
 # 
 # "clus_findCriterionDir <- function(critname)" --
 # 
 # Returns a list with the full criterion name and the subdirectory
 # containing the corresponding test file. Raises an error in case of
 # unknown or ambiguous criterion name.
 # 
 # ------------------------------------------------------------------------
 ##
clus_findCriterionKind <- function(critname) 
{
	# Look among internal criteria
	names <- tolower(getCriteriaNames(TRUE))
	idx <- charmatch(critname, names)
	res <- list()
	if (is.na(idx)) {
		# Look among external criteria
		names <- tolower(getCriteriaNames(FALSE))
		idx <- charmatch(critname, names)
			if (is.na(idx)) {
				stop("\n\t\t>> unknown criterion ",critname)
			} else if (idx == 0) {
				stop("\n\t\t>> ambiguous criterion name ",critname)
			} else {
				res[[1]] <- names[idx]
				res[[2]] <- "extCriteria"
			}
	} else if (idx == 0) {
		stop("\n\t\t>> ambiguous criterion name ",critname)
	} else {
		res[[1]] <- names[idx]
		res[[2]] <- "intCriteria"
	}
	return(res)
}


pkg <- "clusterCrit"

if (require("RUnit", quietly = TRUE)) {
    library(package=pkg, character.only = TRUE)
	
	# Find the criterion name
	res <- clus_findCriterionKind(testName)
	critName <- res[[1]]
	subDir <- res[[2]]
	
	unitTestDir <- system.file("unitTests", package = pkg)
	clusTestDir <- file.path(unitTestDir,subDir)
	testFileName <- paste("runit_",critName,".R",sep="")
	
	# Run the test file
	tst <- runTestFile(file.path(clusTestDir,testFileName), 
		testFuncRegexp = "^test\\..+", 
		verbose = TRUE)

	# Set the output directory
	if (file.access(clusTestDir, 02) != 0) {
		# No write permission, create temporary directory
		tdir <- tempfile(paste(pkg, "unitTests", sep="_"))
		reportsDir <- file.path(tdir, "reports")
	} else {
		reportsDir <- file.path(unitTestDir, "reports")
	}
	dir.create(reportsDir)
	cat("RUnit reports are written into", reportsDir,"\n")

	# Print the results
	printTextProtocol(tst, showDetails = FALSE)
	# Text version
	printTextProtocol(tst, showDetails = TRUE,
					  fileName = file.path(reportsDir,paste("report_", critName, ".txt", sep = "")))
	# HTML version
	printHTMLProtocol(tst, fileName = file.path(reportsDir,paste("report_", critName, ".txt", sep = "")))

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


