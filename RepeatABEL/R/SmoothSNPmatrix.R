#' @title Imputes column means to missing genotypes
#'
#' @description
#' Imputes column means to missing genotypes.
#'
#' @param SNP A matrix including SNP coding.
#' 
#' @author Lars Ronnegard
#' 
SmoothSNPmatrix <-
function(SNP) {
	#Imputes column means to missing values
	m <- ncol(SNP)
    miss <- is.na(colMeans(SNP))
	for (i in (1:m)[miss]){
	SNP[is.na(SNP[,i]), i] <- mean( SNP[!is.na(SNP[ , i]), i] )
	}
	return(SNP)
}
