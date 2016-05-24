

make_random_matrix <- function(range_dim1 = c(50,500), range_dim2 = c(50,500), range_data = c(-1e16,1e16), type="double")
#make_random_matrix <- function(range_dim1 = c(200,1000), range_dim2 = c(200,1000), range_data = c(-1e16,1e16), type="double")
#make_random_matrix <- function(range_dim1 = c(500,1919), range_dim2 = c(1000,5000), range_data = c(-1e16,1e16), type="double")
{
	dim1 <- round(runif(1,range_dim1[1],range_dim1[2]))
	dim2 <- round(runif(1,range_dim2[1],range_dim2[2]))
	data <- runif(dim1*dim2,min=range_data[1],max=range_data[2])
	data <- as(data,type)
	data <- matrix(data,nrow=dim1,ncol=dim2)
	namesCol <- paste("col",c(1:dim2),sep="_")
	namesRow <- paste("row",c(1:dim1),sep="_")
	dimnames(data) <- list(namesRow,namesCol)
	return(data)
}

make_matrix_of_type_with_NA <- function(TYPE="DOUBLE",na.proportion=0.1, range_data, ...)
{
	if (TYPE == "DOUBLE") { 
		if (missing(range_data)) {range_data <- c(-1e100,1e100)} 
		v <- make_random_matrix(range_data=range_data, type = "double", ...);
	} else if (TYPE == "FLOAT") {
		if (missing(range_data)) {range_data <- c(-1e16,1e16)} 
		v <- make_random_matrix(range_data=range_data, type = "double", ...);
	} else if (TYPE == "INT") {
		if (missing(range_data)) {range_data <- c(-2147483648,2147483647-1)} 
		v <- make_random_matrix(range_data=range_data, type = "integer", ...);
	} else if (TYPE == "UNSIGNED_INT") {
		if (missing(range_data)) {range_data <- c(0,4294967295-1)} 
		v <- make_random_matrix(range_data=range_data, type = "integer", ...);
	} else if (TYPE == "SHORT_INT") {
		if (missing(range_data)) {range_data <- c(-32768,32767-1)} 
		v <- make_random_matrix(range_data=range_data, type = "integer", ...);
	} else if (TYPE == "UNSIGNED_SHORT_INT") {
		if (missing(range_data)) {range_data <- c(0,65535-1)} 
		v <- make_random_matrix(range_data=range_data, type = "integer", ...);
	} else if (TYPE == "CHAR") {
		if (missing(range_data)) {range_data <- c(-128,127-1)} 
		v <- make_random_matrix(range_data=range_data, type = "integer", ...);
	} else if (TYPE == "UNSIGNED_CHAR") {
		if (missing(range_data)) {range_data <- c(0,255-1)} 
		v <- make_random_matrix(range_data=range_data, type = "integer", ...);
	} else {
		stop("Unknown TYPE")
	} 
	
	whichNA <- sample(1:length(v),round(length(v)*na.proportion));
	v[whichNA] <- NA;
	
	v
}


checkNumEq <- function(testmatr,test_fv,tolmult=5,quiet=TRUE)
{
	if (!quiet) print("CheckNumEq()");
    matr <- as(test_fv,"matrix")
	if (!quiet) print("testmatr = ");
	if (!quiet) show(testmatr);
	if (!quiet) print("matr = ");
	if (!quiet) show(matr)
	if (!quiet) print("test_fv = ");    
	if (!quiet) show(test_fv)
	checkEqualsNumeric(testmatr,matr,tolerance=tolmult*sqrt(.Machine$double.eps))
	# not 1x1 matrix ?
	if (!is.null(dim(testmatr))) {
	    checkIdentical(dim(testmatr),dim(test_fv))
	}
	dmn <- dimnames(test_fv)
	if (is.null(dmn)) dmn <- get_dimnames(test_fv)
	if (!quiet) cat("dmn=")
	if (!quiet) show(dmn)
	if (!quiet) cat("dimnames(testmatr)=")
	if (!quiet) show(dimnames(testmatr))
	checkIdentical(dimnames(testmatr),dmn)
}	
