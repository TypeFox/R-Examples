# ===========================================================================
# File: "runit_gdi22.R"
#                        Created: 2012-11-13 11:28:57
#              Last modification: 2013-04-20 11:07:09
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@u-paris10.fr>
# Unit test file for the R package clusterCrit.
# ===========================================================================



test.gdi22 <- function() {
	dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsInternal_400_4.Rdata")
	load(file=dataPath, envir=.GlobalEnv)
	idx <- intCriteria(traj_400_4, part_400_4[[4]], c("GDI22"))
	cat(paste("\nFound idx =",idx))
	val <- 15.8152471473492
	cat(paste("\nShould be =",val,"\n"))
	checkEqualsNumeric(idx[[1]],val)
}


