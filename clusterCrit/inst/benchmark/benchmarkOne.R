# ===========================================================================
# File: "runAllTests.R"
#                        Created: 2012-11-01 18:17:51
#              Last modification: 2012-11-13 14:38:33
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
 # "clus_findCritName <- function(critname)" --
 # 
 # Returns a list with the full criterion name and the kind of criterion
 # (internal or external). Raises an error in case of unknown or ambiguous
 # criterion name.
 # 
 # ------------------------------------------------------------------------
 ##
clus_findCritName <- function(critname) 
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
				res[[2]] <- "external"
			}
	} else if (idx == 0) {
		stop("\n\t\t>> ambiguous criterion name ",critname)
	} else {
		res[[1]] <- names[idx]
		res[[2]] <- "internal"
	}
	return(res)
}


pkg <- "clusterCrit"

if (require("rbenchmark", quietly = TRUE)) {
    library(package=pkg, character.only = TRUE)
    
    # Find the criterion name
    res <- clus_findCritName(testName)
    critName <- res[[1]]
    kind <- res[[2]]
    
    if (kind=="internal") {
	dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsInternal_400_4.Rdata")
	load(file=dataPath, envir=.GlobalEnv)
	df <- benchmark(
	    intCriteria(traj_400_4, part_400_4[[4]], c(critName)),
	    columns = c("elapsed", "replications", "relative", "user.self", "sys.self"),
	    replications=100
	)	    
	print(df)
    } else {
    dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsExternal100.Rdata")
    load(file=dataPath, envir=.GlobalEnv)
	df <- benchmark(
	    extCriteria(clus_p2, clus_p3, c(critName)),
	    columns = c("elapsed", "replications", "relative", "user.self", "sys.self"),
	    replications=100
	)	    
	print(df)
    }
    
} else {
    cat("R package 'rbenchmark' cannot be loaded.\n")
}


