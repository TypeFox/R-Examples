# ===========================================================================
# File: "runit_russel_rao.R"
#                        Created: 2012-11-06 20:02:30
#              Last modification: 2015-08-31 10:09:41
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================



test.russel_rao <- function() {
	dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsExternal100.Rdata")
	load(file=dataPath, envir=.GlobalEnv)
	idx <- extCriteria(clus_p2, clus_p3, c("Russel_Rao"))
	val <- 0.170707076787949
	cat(paste("\nShould be =",val,"\n"))
	cat(paste("\nFound idx =",idx,"\n"))
	checkEqualsNumeric(idx[[1]],val)
}


