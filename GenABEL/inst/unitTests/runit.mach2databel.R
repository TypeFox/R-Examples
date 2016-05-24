### --- Test setup ---

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library("RUnit")
    library("GenABEL")
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source("../inst/unitTests/shared_functions.R")
source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.mach2databel <- function()
{
    
#    library("RUnit")
#    library("GenABEL")

    unlink("test*.fv?")

    mli <- read.table("test.mlinfo",head=TRUE,strings=FALSE)
    snpnames <- mli$SNP
    snpnames2 <- rep(0,length(snpnames)*2)
    snpnames2[c(T,F)] <- paste(snpnames,"_11",sep="")
    snpnames2[c(F,T)] <- paste(snpnames,"_01",sep="")
    
    txtdose <- read.table("test.mldose",head=FALSE,strings=FALSE)
    txtdose <- txtdose[,-2]
    txtdose[,1] <- sub("[0-9]*->","",txtdose[,1])
    idnames <- txtdose[,1]
    txtdose <- txtdose[,-1]
    rownames(txtdose) <- idnames
    colnames(txtdose) <- snpnames
    txtdose <- as.matrix(txtdose)

    txtprob <- read.table("test.mlprob",head=FALSE,strings=FALSE)
    txtprob <- txtprob[,-2]
    txtprob[,1] <- sub("[0-9]*->","",txtprob[,1])
    idnames <- txtprob[,1]
    txtprob <- txtprob[,-1]
    rownames(txtprob) <- idnames
    colnames(txtprob) <- snpnames2
    txtprob <- as.matrix(txtprob)
    
    
    fvdose <- mach2databel(imputedg="test.mldose",mlinfo="test.mlinfo",out="test.dose")
    checkIdentical(dim(txtdose),dim(fvdose))
    print(fvdose)
    print(dim(txtdose))
    print(txtdose[1:5,])
    
    
    fvprob <- mach2databel(imputedg="test.mlprob",mlinfo="test.mlinfo",out="test.prob",isprob=TRUE)
	checkIdentical(dim(txtprob),dim(fvprob))
	print(fvprob)
    print(dim(txtprob))
    print(txtprob[1:5,])
    
    checkNumEq(txtdose,fvdose)
    checkNumEq(txtprob,fvprob)
	
	rm(list=ls());gc()
    
    unlink("test*.fv?")    
    
}

