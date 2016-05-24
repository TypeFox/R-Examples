### --- Test setup ---

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library("RUnit")
	library("DatABEL")
}

test.empty <- function(){}
### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

source(paste(path,"/shared_functions.R",sep=""))
source(paste(path,"/mytest_iterator.R",sep=""))
quiet <- TRUE 
	
### --- Test functions ---

test.iterator <- function()
{
#	testmatr <- make_random_matrix()
	
	unlink("tmp*")
	
    testmatr <- matrix(rnorm(10), 2, 5)
    dfo <- as(testmatr, "databel")
	
	# Test sum
    res0 <- apply(testmatr, FUN="sum", MAR=2)
    res1 <- mytest_iterator(dfo, FUN="sum", OUT="R", MAR=2)
    checkEquals(res0, res1)
	
	# Test product
    res2 <- apply(testmatr, FUN="prod", MAR=2)
    res3 <- mytest_iterator(dfo, FUN="prod", OUT="R", MAR=2)
    checkEquals(res2, res3)
	
	# Test sum of powers
	exponent <- 3
	testmatrpow <- testmatr^exponent
	res4 <- apply(testmatrpow, FUN="sum", MAR=2)
	res5 <- mytest_iterator(dfo, FUN="sumpower", OUT="R", MAR=2, POW=exponent)
	checkEquals(res4, res5)
	
	# Test row-wise sum
	res6 <- apply(testmatr, FUN="sum", MAR=1)
	res7 <- mytest_iterator(dfo, FUN="sum", OUT="R", MAR=1)
	checkEquals(res6, res7)
	
	# Test old data type, TO DO!
	#data(srdta)
	#res8 <- apply(srdta@gtdata@gtps, FUN="sum", OUT="R", MAR=2)
	#res9 <- mytest_iterator(srdta, FUN="sum", OUT="R", MAR=2)
	#checkEquals(res8, res9)
	
	# Test filevector output, TO DO!
	res10 <- mytest_iterator(dfo, FUN="sum", OUT="filevector_test.prob", MAR=2)
	
	rm(res10);gc()
	unlink("filevector_test*")
	
	# Test exceptions
    a <- rnorm(100)
    checkException(mytest_iterator(a, FUN="prod", OUT="R", MAR=2))
	checkException(mytest_iterator(dfo, FUN="prod", OUT="R", MAR=3))
    checkException(mytest_iterator(dfo, FUN="non_existent_function", OUT="R", MAR=2))
	checkException(mytest_iterator(dfo, FUN="sumpower", OUT="R", MAR=2, POW="powpow"))

	rm(list=ls());gc()
	
	unlink("tmp*")
}