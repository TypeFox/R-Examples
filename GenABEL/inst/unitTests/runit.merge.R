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

test.merge.bug1641 <- function()
{
	
  require(GenABEL.data)
	data(srdta)
	x1 <- srdta[1:4,1]
	x2 <- srdta[5:10, 2]
	xy <- merge(x1, x2, intersected_snps_only=FALSE)
}

test.merge.bug1676 <- function() 
{
  require(GenABEL.data)
	data(srdta)
	x1 <- srdta[1:4,1]
	x2 <- srdta[5:10, 2]
	xy <- merge(x1,x2)
	phdata(x1)
	phdata(x2)
	phdata(xy)
	checkIdentical(phdata(xy)[1:dim(phdata(x1))[1],],phdata(x1))
	checkIdentical(phdata(xy)[(dim(phdata(x1))[1]+1):(dim(phdata(x1))[1]+dim(phdata(x2))[1]),],phdata(x2))
}
