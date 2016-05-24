### --- Test setup ---
#
# regression test
#

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library(RUnit)
	library(GenABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
#source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.convert.snp.ped <- function()
{
	pref <- ""
	pedf <- paste(pref,"pedin.18",sep="")
	mapf <- paste(pref,"map.18",sep="")
	convert.snp.ped(pedfile=pedf,mapfile=mapf,out="tmp.raw")
	unlink("tmp.raw")
}