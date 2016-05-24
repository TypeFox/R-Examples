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

### load shared functions

#source("../inst/unitTests/shared_functions.R")
source(paste(path,"/shared_functions.R",sep=""))

quiet <- TRUE

### test functions

test.text2databel <- function() 
{
#    library("RUnit")
#    library("DatABEL")
#    source("shared_functions.R")

    unlink("test_matrix*")
    # create matrix
    # data <- make_random_matrix(range_dim1 = c(5,20), range_dim2 = c(5,20))
    data <- make_matrix_of_type_with_NA()
    NR <- dim(data)[1]
    NC <- dim(data)[2]

# create text files
    cat("\nDIM data:",dim(data),"\n")
	write.table(data,file="test_matrix_dimnames.dat",row.names=TRUE,col.names=TRUE,quote=FALSE)
	write.table(data,file="test_matrix_dimnames_TAB.dat",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
	write.table(data,file="test_matrix_colnames.dat",row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(data,file="test_matrix_rownames.dat",row.names=TRUE,col.names=FALSE,quote=FALSE)
    write.table(data,file="test_matrix_NOnames.dat",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write(colnames(data),file="test_matrix.colnames")
    write(rownames(data),file="test_matrix.rownames")
    
# generate identical data and check
    a <- text2databel(infile="test_matrix_dimnames.dat",outfile="test_matrix_dimnames",R_matrix=TRUE) 
    x <- databel("test_matrix_dimnames")
    tmp <- as(x,"matrix")
	if (!quiet) print(data)
	if (!quiet) print(tmp)
    checkEquals(data,tmp)
    data <- tmp

	# test standard options 
    text2databel(infile="test_matrix_dimnames.dat",outfile="test_matrix_dimnames_T",
            R_matrix=TRUE,transpose=TRUE) 
    x <- databel("test_matrix_dimnames_T")
    checkIdentical(t(data),as(x,"matrix"))
    
	# test with TAB separator
	text2databel(infile="test_matrix_dimnames_TAB.dat",outfile="test_matrix_dimnames_TTAB",
			R_matrix=TRUE,transpose=TRUE) 
	x <- databel("test_matrix_dimnames_TTAB")
	checkIdentical(t(data),as(x,"matrix"))
	
	# convert text two filevector format
    
    text2databel(infile="test_matrix_NOnames.dat",outfile="test_matrix_NOnames",
            colnames="test_matrix.colnames",rownames="test_matrix.rownames") 
    x <- databel("test_matrix_NOnames")
    checkIdentical(data,as(x,"matrix"))
    
    text2databel(infile="test_matrix_NOnames.dat",outfile="test_matrix_NOnames_T",
            colnames="test_matrix.colnames",rownames="test_matrix.rownames",transpose=TRUE) 
    x <- databel("test_matrix_NOnames_T")
    checkIdentical(t(data),as(x,"matrix"))
    
    text2databel(infile="test_matrix_NOnames.dat",outfile="test_matrix_NOnames_NOnames") 
    x <- databel("test_matrix_NOnames_NOnames")
    checkEqualsNumeric(data,as(x,"matrix"))
    
    text2databel(infile="test_matrix_NOnames.dat",outfile="test_matrix_NOnames_NOnames_T",transpose=TRUE) 
    x <- databel("test_matrix_NOnames_NOnames_T")
    checkEqualsNumeric(t(data),as(x,"matrix"))
    
    text2databel(infile="test_matrix_rownames.dat",outfile="test_matrix_rownames",
            rownames=1,colnames="test_matrix.colnames") 
    x <- databel("test_matrix_rownames")
    checkIdentical(data,as(x,"matrix"))
    
    text2databel(infile="test_matrix_rownames.dat",outfile="test_matrix_rownames_T",
            rownames=1,colnames="test_matrix.colnames",transpose=TRUE) 
    x <- databel("test_matrix_rownames_T")
    checkIdentical(t(data),as(x,"matrix"))
    
    text2databel(infile="test_matrix_rownames.dat",outfile="test_matrix_rownames_NoColNames",
            rownames=1) 
    x <- databel("test_matrix_rownames_NoColNames")
    checkEqualsNumeric(data,as(x,"matrix"))
    
    text2databel(infile="test_matrix_rownames.dat",outfile="test_matrix_rownames_NoColNames_T",
            rownames=1,transpose=TRUE) 
    x <- databel("test_matrix_rownames_NoColNames_T")
    checkEqualsNumeric(t(data),as(x,"matrix"))
    
    text2databel(infile="test_matrix_colnames.dat",outfile="test_matrix_colnames",
            colnames=1,rownames="test_matrix.rownames") 
    x <- databel("test_matrix_colnames")
    checkIdentical(data,as(x,"matrix"))
    
    text2databel(infile="test_matrix_colnames.dat",outfile="test_matrix_colnames_T",
            colnames=1,rownames="test_matrix.rownames",transpose=TRUE) 
    x <- databel("test_matrix_colnames_T")
    checkIdentical(t(data),as(x,"matrix"))
    
    text2databel(infile="test_matrix_colnames.dat",outfile="test_matrix_colnames_NoRowNames",
            colnames=1) 
    x <- databel("test_matrix_colnames_NoRowNames")
    checkEqualsNumeric(data,as(x,"matrix"))
    
    text2databel(infile="test_matrix_colnames.dat",outfile="test_matrix_colnames_NoRowNames_T",
            colnames=1,transpose=TRUE) 
    x <- databel("test_matrix_colnames_NoRowNames_T")
    checkEqualsNumeric(t(data),as(x,"matrix"))
    
    text2databel(infile="test_matrix_dimnames.dat",outfile="test_matrix_dimnames",R_matrix=TRUE) 
    x <- databel("test_matrix_dimnames")
    checkIdentical(data,as(x,"matrix"))
    
    text2databel(infile="test_matrix_dimnames.dat",outfile="test_matrix_dimnames_T",
            R_matrix=TRUE,transpose=TRUE) 
    x <- databel("test_matrix_dimnames_T")
    checkIdentical(t(data),as(x,"matrix"))
    
    text2databel(infile="test_matrix_NOnames.dat",outfile="test_matrix_sub",
            skipcols=3,skiprows=2) 
    x <- databel("test_matrix_sub")
    checkEqualsNumeric(data[c(3:dim(data)[1]),c(4:dim(data)[2])],as(x,"matrix"))
    
    # stupid extended matrix in non-R format
    newmat <- make_matrix_of_type_with_NA(range_dim1 = c(NR+2,NR+2), range_dim2 = c(NC+3,NC+3))
    newmat[3:(NR+2),4:(NC+3)] <- data
    newmat[2,4:(NC+3)] <- colnames(data)
    newmat[3:(NR+2),3] <- rownames(data)
    newmat
    write.table(newmat,file="test_matrix_strange.dat",col.names=FALSE,row.names=FALSE,quote=FALSE)
    
    text2databel(infile="test_matrix_strange.dat",outfile="test_matrix_strange",
            colnames=2,rownames=3) 
    x <- databel("test_matrix_strange")
    checkEquals(data,as(x,"matrix"))

	rm(list=ls());gc()
	
	unlink("test_matrix*")
    
}