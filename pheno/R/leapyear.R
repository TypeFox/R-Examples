# TRUE if leap year, FALSE otherwise
leapyear <- function(y) {
	
	res <- .C("Cleapyear",year=as.integer(y),ly=integer(1),PACKAGE="pheno")

	if(res$ly==1) { return(TRUE)}
	else { return(FALSE)}
}
