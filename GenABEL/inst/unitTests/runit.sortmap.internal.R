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

test.sortmap.internal.bug1673 <- function()
{
	nmrkInChr <- c(7,5,3)
	chrLevels <- c("cH1","X","Y")
	chroms <- mappos <- rep(NA,sum(nmrkInChr))
	k <- 1
	for (i in 1:length(nmrkInChr)) {
		for (j in seq(nmrkInChr[i],1,-1)) {
			chroms[k] <- chrLevels[i]
			mappos[k] <- j
			k <- k + 1
		}
	}
	a <- sortmap.internal(chroms,mappos)
	checkTrue( all(is.finite(a$chnum)) )
}
