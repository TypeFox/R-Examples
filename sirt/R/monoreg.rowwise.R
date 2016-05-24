
##############################################################
# monotone regression for all rows in a matrix
monoreg.rowwise <- function(yM,wM){ 
    yM <- as.matrix(yM)
	wM <- as.matrix(wM)
	res <- .Call("monoreg_rowwise_Cpp", yM , wM , PACKAGE = "sirt")
	return(res)
					}	
##############################################################
# monotone regression for all columns in a matrix
monoreg.colwise <- function(yM,wM){
    yM <- as.matrix(yM)
	wM <- as.matrix(wM)
	res <- .Call("monoreg_rowwise_Cpp", t(yM) , t(wM) , PACKAGE = "sirt")
	return(t(res))
					}	