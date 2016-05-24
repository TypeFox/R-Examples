# calculates maximal daylength maxdl [h] at a certain latitude lat [degrees]
maxdaylength <- function(l) {

	res <- .C("Cmaxdaylength",l=as.double(l),maxdl=double(1),PACKAGE="pheno")

	return(res$maxdl)
}
