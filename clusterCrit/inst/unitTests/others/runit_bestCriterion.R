# ===========================================================================
# File: "runit_bestCriterion.R"
#                        Created: 2012-11-06 20:02:30
#              Last modification: 2012-11-13 11:46:11
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================



test.bestCriterion <- function() {
	
	dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsInternal_400_4.Rdata")
	load(file=dataPath, envir=.GlobalEnv)

	crits <- getCriteriaNames(TRUE)
	len <- length(crits)
	vals <- matrix(nrow=len,ncol=6)
	rownames(vals) <- crits
	colnames(vals) <- paste("P",2:7,sep="")

	for (k in 2:7) {
		traj <- traj_400_4
		part <- part_400_4[[k]]
		
		# Calculate all the internal indices
		indices <- intCriteria(traj,part,"all")
		
		for (i in 1:len) {
			vals[i,k-1] <- as.numeric(indices[[i]])
		}
	}

	result <- vector(mode="numeric", length=len)
	names(result) <- crits
	
	for (i in 1:len) {
		cn <- crits[i]
			idx <- bestCriterion(vals[i,],cn)
			if (!is.null(idx)) {
				result[i] <- idx+1
			}
	}
	cat("\nBest criterion yields (expected 4):\n")
	cat(result,"\n")
	
	expected <- c(3,7,4,4,5,4,4,2,7,4,4,4,4,4,4,4,4,4,4,4,4,2,2,2,4,4,4,5,4,2,4,3,7,7,3,4,4,7,4,3,4,4)
	checkEqualsNumeric(result,expected)
}


