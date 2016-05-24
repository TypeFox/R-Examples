
########################################################
# Call to Rcpp function
ma.scale2 <- function (x , missings =FALSE ){ 
    x_ <- as.matrix(x)
	if ( ! missings ){
		res <- .Call("scale2_C", x_, PACKAGE = "miceadds")
					} else {
		res <- .Call("scale2_NA_C", x_, PACKAGE = "miceadds")
					}
	colnames(res) <- colnames(x)
	return(res)
					}
##########################################################