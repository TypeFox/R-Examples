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

test.convert.snp <- function()
{
	# library(GenABEL)
	# check equivalence of convertion
	dta <- read.table("gatest.tped",sep="\t",head=FALSE,strings=FALSE)
	dta <- dta[,-c(1:4)]
	dta <- t(dta)
	dta <- sub(" ","/",dta)
	dta <- sub("0/0",NA,dta)
	dimnames(dta) <- NULL

	convert.snp.tped(tped="gatest.tped", tfam="gatest.tfam", out="gatest.raw")
	data <- load.gwaa.data(gen="gatest.raw", phe="gatest.phe")
	dta1 <- as.character(data)
	dimnames(dta1) <- NULL
	checkEquals(dta,dta1)
	
	convert.snp.ped(pedfile="gatest.ped", mapfile="gatest.map", out="gatest.raw",mapHasHeaderLine=FALSE)
	data <- load.gwaa.data(gen="gatest.raw", phe="gatest.phe")
	dta1 <- as.character(data)
	dimnames(dta1) <- NULL
	checkEquals(dta,dta1)
	
	convert.snp.illumina(infile="gatest.illu", out="gatest.raw")
	data <- load.gwaa.data(gen="gatest.raw", phe="gatest.phe")
	dta1 <- as.character(data)
	dimnames(dta1) <- NULL
	checkEquals(dta,dta1)
	
	# check possibility of merge
	
	convert.snp.tped(tped="gatest1.tped", tfam="gatest1.tfam", out="gatest1.raw", strand="+")
	data1 <- load.gwaa.data(gen="gatest1.raw", phe="gatest1.phe")
	data1 <- data1[,sample(1:nsnps(data1),nsnps(data1))]
	coding(data1)
	mrg <- merge(gtdata(data),gtdata(data1))
	mrg1 <- merge(data,data1)
	checkEquals(as.character(mrg$data),as.character(gtdata(mrg1)))
	
	# here evaluate merged SNP one by one?

	unlink("*.raw")
}