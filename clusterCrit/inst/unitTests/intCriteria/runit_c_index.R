# ===========================================================================
# File: "runit_c_index.R"
#                        Created: 2012-11-13 11:28:57
#              Last modification: 2015-08-31 10:09:41
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@u-paris10.fr>
# Unit test file for the R package clusterCrit.
# ===========================================================================



test.c_index <- function() {
	dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsInternal_400_4.Rdata")
	load(file=dataPath, envir=.GlobalEnv)
	idx <- intCriteria(traj_400_4, part_400_4[[4]], c("C_index"))
	cat(paste("\nFound idx =",idx))
	val <- 7.0592193043113e-06
	cat(paste("\nShould be =",val,"\n"))
	checkEqualsNumeric(idx[[1]],val)
}


