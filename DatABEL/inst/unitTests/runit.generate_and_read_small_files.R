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

test.as.databel <- function()
{
#    library("RUnit")
#    library("DatABEL")
#    source("shared_functions.R")

    unlink("tmp*")


    testmatr <- make_matrix_of_type_with_NA()
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    disconnect(test_fv)
    connect(test_fv)
    checkNumEq(testmatr,test_fv,quiet=quiet)
    

    test_fv <- matrix2databel(testmatr,file="tmp",type="DOUBLE")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    disconnect(test_fv)
    connect(test_fv)
    checkNumEq(testmatr,test_fv,quiet=quiet)
    

    # try to ini using the same name
    # print(system("ls -al"))
    checkException(test_fv <- matrix2databel(testmatr,file="tmp",type="DOUBLE"))
    
    rm(test_fv);gc()
    unlink("tmp*")

	if (!quiet) print("A4");

    #define UNSIGNED_SHORT_INT 1
    testmatr <- make_matrix_of_type_with_NA("UNSIGNED_SHORT_INT")
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)

    test_fv <- matrix2databel(testmatr,file="tmp",type="UNSIGNED_SHORT_INT")
    checkNumEq(testmatr,test_fv,quiet=quiet)

    rm(test_fv,testmatr);gc()
    unlink("tmp*")

	if (!quiet) print("A5");

    #define SHORT_INT          2
	testmatr <- make_matrix_of_type_with_NA("SHORT_INT")
	test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)

    test_fv <- matrix2databel(testmatr,file="tmp",type="SHORT_INT")
    checkNumEq(testmatr,test_fv,quiet=quiet)

    rm(test_fv,testmatr);
    gc()
    unlink("tmp*")
    

    #define UNSIGNED_INT       3
    testmatr <- make_matrix_of_type_with_NA("UNSIGNED_INT")
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)


    test_fv <- matrix2databel(testmatr,file="tmp",type="UNSIGNED_INT")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    
    rm(test_fv,testmatr);gc()
    unlink("tmp*")

	if (!quiet) print("A7");

    #define INT                4
    testmatr <- make_matrix_of_type_with_NA("INT")
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    
    test_fv <- matrix2databel(testmatr,file="tmp",type="INT")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    
    rm(test_fv,testmatr);gc()
    unlink("tmp*")
    
    #define FLOAT              5
    testmatr <- make_matrix_of_type_with_NA("FLOAT")
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)
	
    test_fv <- matrix2databel(testmatr,file="tmp",type="FLOAT")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    
    rm(test_fv,testmatr);gc()
    unlink("tmp*")
    
    #define DOUBLE             6
    testmatr <- make_matrix_of_type_with_NA("DOUBLE")
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)

    test_fv <- matrix2databel(testmatr,file="tmp",type="DOUBLE")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    rm(test_fv);gc()
    unlink("tmp*")

    #define CHAR              7
    testmatr <- make_matrix_of_type_with_NA("CHAR")
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)

    test_fv <- matrix2databel(testmatr,file="tmp",type="CHAR")
    checkNumEq(testmatr,test_fv,quiet=quiet)
    rm(test_fv);gc()
    unlink("tmp*")

    #define UNSIGNED_CHAR            8
    testmatr <- make_matrix_of_type_with_NA("UNSIGNED_CHAR")
    test_fv <- as(testmatr,"databel")
    checkNumEq(testmatr,test_fv,quiet=quiet)


    rm(test_fv,testmatr);gc()
	rm(list=ls());gc()
    unlink("tmp*")
}

test.save_as <- function()
{

#    library("RUnit")
#    library("DatABEL")
#    source("shared_functions.R")

    unlink("tmp*")
    testmatr <- make_matrix_of_type_with_NA(range_data=c(-1000,1000))
	if (!quiet) print(testmatr)
    dba <- matrix2databel(testmatr,file="tmp")
    dfa <- as(dba,"databel")
	if (!quiet) print(as.matrix(dba))
    checkNumEq(testmatr,dba,quiet=quiet)
	if (!quiet) print(as(dfa,"matrix"))
    checkNumEq(testmatr,dfa,quiet=quiet)
    
    dims <- dim(testmatr)
    subrows <- sample(1:dims[1],2+floor(runif(1,min=(dims[1]-3)/10,max=(dims[1]-3))))
    subcols <- sample(1:dims[2],2+floor(runif(1,min=(dims[2]-3)/10,max=(dims[2]-3))))
    sub_tm <- testmatr[subrows,subcols]
	if (!quiet) print(subrows)
	if (!quiet) print(subcols)
	if (!quiet) print(sub_tm)
    sub_dba <- save_as(dba[subrows,subcols],file="tmp_sub_1")
    checkNumEq(sub_tm,sub_dba,quiet=quiet)
    sub_dba <- save_as(dba,row=subrows,col=subcols,file="tmp_sub_2")
	if (!quiet) print(as.matrix(sub_dba))
	if (!quiet) print(as.matrix(as(sub_dba,"databel_filtered_R")))
    checkNumEq(sub_tm,sub_dba,quiet=quiet)
    
    sub_dfa <- save_as(dfa[subrows,subcols],file="tmp_dfa_sub")
	if (!quiet) print(sub_tm)
	if (!quiet) print(as(dfa[subrows,subcols],"matrix"))
	if (!quiet) print(as(sub_dfa,"matrix"))
    checkNumEq(sub_tm,as(sub_dfa,"matrix"),quiet=quiet)    
    
    sub_dfa <- save_as(dfa,row=subrows,col=subcols,file="tmp_dfa_sub1")
	if (!quiet) print(sub_tm)
	if (!quiet) print(as(dfa[subrows,subcols],"matrix"))
	if (!quiet) print(as(sub_dfa,"matrix"))
    checkNumEq(sub_tm,as(sub_dfa,"matrix"),quiet=quiet)
	
	rm(list=ls());gc()
    
    unlink("tmp*")
}

