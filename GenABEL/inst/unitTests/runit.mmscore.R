### --- Test setup ---
#
# regression tests
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

test.mmscore <- function()
{
  require(GenABEL.data)
	data(ge03d2.clean)
	dta <- ge03d2.clean[1:200,autosomal(ge03d2.clean)[1:3000]]
	gkin <- ibs(dta[,1:2000],w="freq")
	h2 <- polygenic(height,data=dta,kin=gkin)
	xNew <- mmscore(h2,dta)
	xNew_results <- results(xNew)


# The test_mmscore_results_revision_1321 contains results from mmscore of SNV revision 1321. Let's assume that those results are correct and use it the future as a reference.
#	test_mmscore_results_revision_1321 <- results(xNew)
#	write.table(test_mmscore_results_revision_1321, file="../exdata/test_mmscore_results_revision_1321", quote=F)
	
	test_mmscore_results_revision_1321 <- read.table("test_mmscore_results_revision_1321", stringsAsFactors=F)
	
	test_mmscore_results_revision_1321$Chromosome <- as.factor(test_mmscore_results_revision_1321$Chromosome)
	test_mmscore_results_revision_1321$Strand <- as.factor(test_mmscore_results_revision_1321$Strand)
	test_mmscore_results_revision_1321$A1 <- as.factor(test_mmscore_results_revision_1321$A1)
	test_mmscore_results_revision_1321$A2 <- as.factor(test_mmscore_results_revision_1321$A2)
		
	checkEquals(xNew_results, test_mmscore_results_revision_1321)
}
