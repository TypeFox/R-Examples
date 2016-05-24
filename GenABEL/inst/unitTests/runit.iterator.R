### --- Test setup ---

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library(RUnit)
	library(GenABEL)
	library(DatABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.summary_snp_data <- function()
{
  require(GenABEL.data)
	data(ge03d2)
	srdta <- ge03d2 #[,autosomal(ge03d2)]
	
	t0 <- proc.time()
	smrOld <- summary.snp.data_old(gtdata(srdta))
	tOld <- proc.time() - t0
#smrOld
	
	srdtaOLDT <- srdta
#gtdata(srdta)
	t0 <- proc.time()
	smrNewOLDT <- summary(gtdata(srdtaOLDT))
	tNewOLDT <- proc.time() - t0
#smrNewOLDT
	checkIdentical(smrOld,smrNewOLDT)
	
	srdtaMATRIX <- srdta
	srdtaMATRIX@gtdata@gtps <- as.numeric(srdta)
#gtdata(srdtaMATRIX)
	t0 <- proc.time()
	smrNewMATRIX <- summary(gtdata(srdtaMATRIX))
	tNewMATRIX <- proc.time() - t0
#smrNewMATRIX
	checkIdentical(smrOld,smrNewMATRIX)
	
	srdtaDA <- srdta
	srdtaDA@gtdata@gtps <- matrix2databel(as.numeric(gtdata(srdta)),
			file="tmp",cache=16,type="UNSIGNED_CHAR",readonly=TRUE)
#gtdata(srdtaDA)
	t0 <- proc.time()
	smrNewDA <- summary(gtdata(srdtaDA))
	tNewDA <- proc.time() - t0
#smrNewDA
	checkIdentical(smrOld,smrNewDA)
	
	print(tOld)
	print(tNewOLDT)
	print(tNewMATRIX)
	print(tNewDA)
	
	rm(list=ls())
	gc()
	unlink("tmp.fv?")
}
