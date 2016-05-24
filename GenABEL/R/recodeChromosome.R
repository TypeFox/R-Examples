#' Change chromosomal coding
#' 
#' Recoding of chromosomes according to the provided 'rules' for from -> to 
#' pairs. Most common use is anticipated when importing data from other 
#' software using only integers to represent chromosomes. In this situation 
#' the list of rules may look like this: list(24="X",25="Y",26="mt").
#' 
#' @note 'from' entries should be unique and not overlap with entries in 'to'
#' 
#' @param data object of class for which 'chromosome' method is defined, e.g. 
#' 'gwaa.data', 'snp.data', 'scan.gwaa'
#' 
#' @param rules list of pairs 'from=to'; the chromosomes coded in the original data 
#' set with 'from' will be recoded with 'to' value
#' 
#' @param quiet if summary of recoding should not be printed to the screen
#' 
#' @return modified 'data' object
#' 
#' @author Yurii Aulchenko
#' 
#' @examples 
#' data(ge03d2)
#' table(chromosome(ge03d2))
#' # merge chromosome 3 and X, call chromosome 2 as 15
#' newdat <- recodeChromosome(ge03d2,rules=list("3"="X","2"=15))
#' table(chromosome(ge03d2),chromosome(newdat))
#' 
#' @keywords manip
#' 
#' 

recodeChromosome <- function(data,rules,quiet=FALSE) 
{
	chrom <- try( chromosome ( data ) )
	if (class(chrom) == "try-error") stop("'data' should have gwaa.data or snp.data class (1)")
	if (class(rules) != "list") stop("'rules' should be a list")
	if (anyDuplicated(names(rules))) stop("duplicated entries in 'rules' from-entries")
	if (any( names(rules) %in% unlist(rules) )) stop("overlap in 'rules' between from and to entries")
	for (fromName in names(rules)) {
		toName <- rules[[fromName]]
		if (length(toName) != 1) 
			stop(paste('rules list element with name',fromName,'has #entries <> 1'))
		saveOpt <- getOption("warn")
		options("warn" = -1)
		if (!(toName == as.integer(toName) | toName %in% c("X","Y","mt") ) ) {
			warning(paste("to-name",toName,"is neither integer nor one of 'X', 'Y', 'mt'"))
		}
		options("warn" = saveOpt)
		toBeRecoded <- which(chrom == as.character(fromName))
		chrom[toBeRecoded] <- as.character(toName)
		if (length(toBeRecoded) >= 1) {
			cat("Recoded chromosome for",length(toBeRecoded),"SNPs (",fromName,"->",toName,")\n")
		} else if (!quiet) {
			cat("No chromosome coded as",fromName,"found\n")
		}
	}
	if (class(data) == "gwaa.data") 
		data@gtdata@chromosome <- as.factor(chrom)
	else if (class(data) == "snp.data")
		data@chromosome <- as.factor(chrom)
	else if (class(data) == "scan.gwaa")
		data@annotation[,"Chromosome"] <- as.factor(chrom)
	else 
		stop("'data' should have gwaa.data or snp.data class (2)")
	return(data)
}