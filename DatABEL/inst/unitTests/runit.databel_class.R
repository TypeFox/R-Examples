### --- Test setup ---

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library("RUnit")
    library("DatABEL")
}

#test.empty <- function(){}
### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source("../inst/unitTests/shared_functions.R")
source(paste(path,"/shared_functions.R",sep=""))

quiet <- TRUE 

### --- Test functions ---

#
# add NAs to the test
#

test.databel_class <- function()
{
    cat ("Running test.databel_class\n");
#    library("RUnit")
#    library("DatABEL")
#    source("shared_functions.R")
    
    unlink("tmp*.fv?")
    unlink("test*.fv?")
    
	# would fail with full range !!!
#	dta <- make_random_matrix(range_data=c(-10,10))
	dta <- make_matrix_of_type_with_NA()
	null <- as(dta,"databel")
    checkNumEq(dta,null,quiet=quiet)
    checkIdentical("databel",class(null)[1])
    x <- null
    checkIdentical("databel",class(x)[1])
    checkNumEq(dta,x,quiet=quiet)
    x1 <- x
    checkIdentical(as(x,"matrix"),as(null,"matrix"))
    checkNumEq(dta,x1,quiet=quiet)
    disconnect(x)
	checkNumEq(dta,null,quiet=quiet)
	checkNumEq(dta,x1,quiet=quiet)
	connect(x)
    checkIdentical("databel",class(x)[1])
    checkNumEq(dta,x,quiet=quiet)
    rm(null);gc()
    checkNumEq(dta,x,quiet=quiet)
    rm(x1);gc()
    x1 <- x[]
    disconnect(x)
    checkNumEq(dta,x1,quiet=quiet)
    connect(x)
    rm(x1);gc()
    x1 <- x
    disconnect(x)
    checkNumEq(dta,x1,quiet=quiet)
    connect(x)
    rm(x1);gc()
    x1 <- x
    rm(x);gc()
    checkNumEq(dta,x1,quiet=quiet)
    rm(x1);gc()
    
    x <- as(dta,"databel")
    checkIdentical("databel",class(x)[1])
    checkNumEq(dta,x,quiet=quiet)
    disconnect(x)
    connect(x)
    checkIdentical("databel",class(x)[1])
    checkNumEq(dta,x,quiet=quiet)
    x1 <- x
    checkNumEq(dta,x1,quiet=quiet)
    
    
    # make identical
#    dta <- as(x,"matrix")
#    testobject(dta,as(x,"matrix"))
    
# submatrix
    rCol <- sample(1:ncol(x),ceiling(ncol(x)/10))
    rCol
    rRow <- sample(1:nrow(x),ceiling(nrow(x)/10))
    rRow
    
    dta1 <- dta[rRow,rCol,drop=FALSE]
    x1 <- x[rRow,rCol]
    checkNumEq(dta1,x1,quiet=quiet)
    disconnect(x1)
    checkNumEq(dta1,x1,quiet=quiet)
    connect(x1)
	
    dta1 <- dta[sort(rRow),sort(rCol),drop=FALSE]
    x1 <- x[sort(rRow),sort(rCol)]
    checkNumEq(dta1,x1,quiet=quiet)
    
    rCol <- sample(1:ncol(x1),ceiling(ncol(x)/10))
    rCol
    rRow <- sample(1:nrow(x1),ceiling(nrow(x)/10))
    rRow
    
    dta2 <- dta1[rRow,rCol,drop=FALSE]
    x2 <- x1[rRow,rCol]
    checkNumEq(dta2,x2,quiet=quiet)
    
    logCol <- (runif(ncol(dta))<=0.1)
    while (!any(logCol)) logCol <- (runif(ncol(dta))<=0.1)
    logRow <- (runif(nrow(dta))<=0.1)
    while (!any(logRow)) logRow <- (runif(nrow(dta))<=0.1)
    table(logCol)
    table(logRow)
    
    dta1 <- dta[logRow,logCol,drop=FALSE]
    x1 <- x[logRow,logCol,drop=FALSE]
### would fail with checkIdentical
#    checkIdentical(dta1,as(x1,"matrix"))
    checkNumEq(dta1,x1,quiet=quiet)
    
    logCol <- (runif(ncol(dta1))<=0.1)
    while (!any(logCol)) logCol <- (runif(ncol(dta1))<=0.1)
    logRow <- (runif(nrow(dta1))<=0.1)
    while (!any(logRow)) logRow <- (runif(nrow(dta1))<=0.1)
    table(logCol)
    table(logRow)

    dta2 <- dta1[logRow,logCol,drop=FALSE]
    x2 <- x1[logRow,logCol]
### would fail with checkIdentical
#	checkIdentical(dta2,as(x2,"matrix"))
    checkNumEq(dta2,x2,quiet=quiet)

    dta1 <- dta[c(T,F,F,F,F),c(F,F,F,F,T),drop=FALSE]
    dim(x)
    x1 <- x[c(T,F,F,F,F),c(F,F,F,F,T)]
### would fail with checkIdentical
#    checkIdentical(dta1,as(x1,"matrix"))
    checkNumEq(dta1,x1,quiet=quiet)

# submatrix and permute
# on integer 
# on logical
# on names
    
    apply(dta1,FUN=sum,MAR=1)
#apply(x1,FUN=sum,MAR=1)
    
    dn <- dimnames(x)
    ndn <- dn
    ndn[[1]] <- sample(dn[[1]],ceiling(length(dn[[1]])/3))
    ndn[[2]] <- sample(dn[[2]],ceiling(length(dn[[2]])/5))
    ndn
    dta1 <- dta[ndn[[1]],ndn[[2]],drop=FALSE]
    x1 <- x[ndn[[1]],ndn[[2]]]
    dta1
    x1
### would fail with checkIdentical
#    checkIdentical(dta1,as(x1,"matrix"))
    checkNumEq(dta1,x1,quiet=quiet)
	
	rm(list=ls());gc()
    
    unlink("tmp*.fv?")
    unlink("test*.fv?")
    
}

test.databel_class_setReadOnly <- function()
{
	cat ("Running test.databel_class_setReadOnly\n");    
	
	unlink("tmp*.fv?")
	unlink("test*.fv?")
	
	dta <- make_matrix_of_type_with_NA()
	dfo <- as(dta,"databel")
	checkNumEq(dta,dfo,quiet=quiet)
	
	dfo[1,1] <- 11
	checkEquals(11,as(dfo[1,1],"vector"))
	setReadOnly(dfo) <- TRUE
	checkException(dfo[1,1] <- 12)
	checkEqualsNumeric(11,as(dfo[1,1],"vector"))
	setReadOnly(dfo) <- FALSE
	dfo[1,1] <- 12
	checkEqualsNumeric(12,as(dfo[1,1],"vector"))
	

	rm(list=ls());gc()
	
	unlink("tmp*.fv?")
	unlink("test*.fv?")
	
}