#'
#' Access to slots of \link{gwaa.data-class} object
#' 
#' Functions to access slots of \link{gwaa.data-class} object. In particular:
#' 
#' phdata: phenotypic data
#' 
#' gtdata: genotypic data
#' 
#' chromosome: SNP chromosome
#' 
#' coding: SNP coding
#' 
#' map: SNP map position
#' 
#' strand: SNP strand
#' 
#' male: subject's sex (who is male?)
#' 
#' idnames: subject IDs
#' 
#' snpnames: marker IDs
#' 
#' nids: number of subjects
#' 
#' nsnps: number of markers
#' 
#' @usage phdata(data)
#' gtdata(data) 
#' chromosome(data)
#' coding(data)
#' idnames(data)
#' male(data)
#' map(data)
#' nids(data)
#' nsnps(data)
#' snpnames(data)
#' strand(data)
#' 
#' @aliases gtdata chromosome coding gtdata idnames male map nids nsnps phdata snpnames strand
#' 
#' @param data object of \link{gwaa.data-class}
#' 
#' @return content of the particular slot
#'

phdata <- function(data) {
	if (class(data) == "gwaa.data") return(data@phdata)
	else stop("data should be of class 'gwaa.data'")
}