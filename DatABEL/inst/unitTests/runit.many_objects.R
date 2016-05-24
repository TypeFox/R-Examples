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

test.many_objects <- function()
{
    unlink("tmp*.fv?")
    unlink("test*.fv?")
    
    madta <- make_matrix_of_type_with_NA(range_dim1=c(25,100),range_dim2=c(25,100))
    dadta <- matrix2databel(madta,file="tmp",cache=1)
    dal <- list()
    for (i in c(1:256)) {
        dal[[i]] <- dadta[c(1:(dim(madta)[1]-1)),c(1:(dim(madta)[2]-1))]
        checkNumEq(madta[c(1:(dim(madta)[1]-1)),c(1:(dim(madta)[2]-1))],dal[[i]],quiet=quiet)
        checkNumEq(madta[c(1:(dim(madta)[1]-1)),c(1:(dim(madta)[2]-1))],
                dadta[c(1:(dim(madta)[1]-1)),c(1:(dim(madta)[2]-1))],quiet=quiet)
        checkEquals(madta[c(1:(dim(madta)[1]-1)),c(1:(dim(madta)[2]-1))],
                as(dal[[i]],"matrix"))
    }
    rm(dal,dadta);gc()
    
    unlink("tmp*.fv?")
    unlink("test*.fv?")
    
}