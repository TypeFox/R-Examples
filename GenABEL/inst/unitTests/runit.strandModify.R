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

test.strandModify <- function()
{
  require(GenABEL.data)
	data(ge03d2ex)
	str <- strand(ge03d2ex)
	dta1 <- ge03d2ex
	uStrand <- rep("u",nsnps(dta1))
	strand(dta1) <- uStrand
	dta2 <- patch_strand(data=dta1,snpid=snpnames(dta1),strand=uStrand)
	table(strand(dta1),strand(ge03d2ex))
	table(strand(dta2),strand(ge03d2ex))
	strand(dta1) <- str
	dta2 <- patch_strand(data=dta2,snpid=snpnames(dta2),strand=str)
	checkIdentical(dta1,ge03d2ex)
	checkIdentical(dta2,ge03d2ex)
	table(strand(dta1),strand(ge03d2ex))
	table(strand(dta2),strand(ge03d2ex))
}
