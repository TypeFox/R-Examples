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

test.recodeChromosome <- function()
{
  require(GenABEL.data)
	data(ge03d2ex)
	lst <- list("3"="X","2"=15)
	dta0 <- chromosome(recodeChromosome(ge03d2ex,lst))
	dta1 <- chromosome(recodeChromosome(gtdata(ge03d2ex),lst))
	qts <- qtscore(dm2,ge03d2ex)
	dta2 <- chromosome(recodeChromosome(qts,lst))
	checkIdentical(dta0,dta1)
	checkIdentical(dta0,dta2)
	
	checkException( dta <- recodeChromosome(10,lst) )
	checkException( dta <- recodeChromosome(ge03d2ex,c(10,12)) )
	lst <- list("1"="X","1"="Y")
	checkException( dta <- recodeChromosome(ge03d2ex,lst) )
	lst <- list("1"="2","2"="3")
	checkException( dta <- recodeChromosome(ge03d2ex,lst) )
	lst <- list("1"="oppa","2"=c(1,2))
	checkException( dta <- recodeChromosome(ge03d2ex,lst) )
	lst <- list("1"="0","2"=c())
	checkException( dta <- recodeChromosome(ge03d2ex,lst) )
}
