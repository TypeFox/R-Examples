# ===========================================================================
# File: "runit_rogers_tanimoto.R"
#                        Created: 2012-11-06 20:02:30
#              Last modification: 2015-08-31 09:58:29
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================



test.rogers_tanimoto <- function() {
	dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsExternal100.Rdata")
	load(file=dataPath, envir=.GlobalEnv)
	idx <- extCriteria(clus_p2, clus_p3, c("Rogers_Tanimoto"))
	val <- 0.344195519348269
	cat(paste("\nShould be =",val,"\n"))
	cat(paste("\nFound idx =",idx,"\n","\n"))
	checkEqualsNumeric(idx[[1]],val)
}


