

make_random_matrix <- function(range_dim1 = c(200,1000), range_dim2 = c(200,1000), range_data = c(-1e16,1e16), type="double")
#make_random_matrix <- function(range_dim1 = c(500,1919), range_dim2 = c(1000,5000), range_data = c(-1e16,1e16), type="double")
#make_random_matrix <- function(range_dim1 = c(1000,3000), range_dim2 = c(1000,3000), range_data = c(-1e16,1e16), type="double")
#make_random_matrix <- function(range_dim1 = c(100,300), range_dim2 = c(100,300), range_data = c(-1e16,1e16), type="double")
{
	dim1 <- round(runif(1,range_dim1[1],range_dim1[2]))
	dim2 <- round(runif(1,range_dim2[1],range_dim2[2]))
	data <- runif(dim1*dim2,range_data[1],range_data[2])
	data <- as(data,type)
	data <- matrix(data,nrow=dim1,ncol=dim2)
	namesCol <- paste("col",c(1:dim2),sep="_")
	namesRow <- paste("row",c(1:dim1),sep="_")
	dimnames(data) <- list(namesRow,namesCol)
	return(data)
}

checkNumEq <- function(testmatr,test_fv,tolmult=5)
{
	checkIdentical(dim(testmatr),dim(test_fv))
	checkEqualsNumeric(testmatr,as(test_fv,"matrix"),tolerance=tolmult*sqrt(.Machine$double.eps))
	dmn <- dimnames(test_fv)
	if (is.null(dmn)) dmn <- get_dimnames(test_fv)
	checkIdentical(dimnames(testmatr),dmn)
}	
