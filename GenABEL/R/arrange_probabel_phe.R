#'
#' arranges ProbABEL phenotype-file
#' 
#' Function to arrange ProbABEL phenotype-file; it takes 
#' phenotypic data as input and 
#' aligns that with genotypic data of ProbABEL
#' 
#' @param modelterms vector of character, which  
#' specifies the variables to be included into 
#' ProbABEL phenotype-file. Should contain, 
#' and start with 'id' column, which should 
#' provide the same ID codes as these in gendata
#' 
#' @param phedata phenotypic data (matrix, data.frame, or 
#' \link{gwaa.data-class} object)
#' 
#' @param gendata genetic data to be used with ProbABEL,
#' either databel-class object or name of the 
#' (index/data) file containing filevector data for 
#' ProbABEL
#' 
#' @param file name of the ProbABEL phenotype file
#' 
#' @return file with phenotypes ready for use with 
#' ProbABEL
#' 
#' @author Yurii Aulchenko
#' 
#' @keywords IO manip
#'

arrange_probabel_phe <- function(modelterms,phedata,gendata,file="probabel.PHE")
{
	if (class(phedata) == "gwaa.phedata") phedata <- phdata(phedata)
	else if (is.character(phedata)) 
		phedata <- read.table(phedata,header=TRUE,stringsAsFactors=FALSE)	
	else if (class(phedata) != "data.frame" && class(phedata) != "matrix")
		stop("phedata should be of class 'gwaa.phedata', 'data.frame', or 'matrix'")
	
	if (missing(modelterms)) 
		stop("you need to specify model terms (vector of names)")
	
	if (!any(colnames(phedata)=="id")) 
		stop("phedata should contain 'id'")
	
	if (any(is.na(match(modelterms,colnames(phedata))))) 
		stop("some modelterms are missing from phedata")
	
	if (is.character(gendata)) { 
		if (require(DatABEL)) 
			gendata <- databel(gendata)
		else 
			stop ("this function requires DatABEL package to be installed")
	} else if (class(gendata) == "databel") gendata <- databel(gendata)
	else if (class(gendata) == "databel") {}
	else stop("gendata should be 'databel-class' object or name of *.FV? file")
	
	dosephe <- data.frame(id=dimnames(gendata)[[1]],stringsAsFactors=F)
	phe <- merge(dosephe,phedata,by="id",all.x=T,all.y=F)
	rownames(phe) <- phe$id
	phe <- phe[as.character(dosephe$id),modelterms]
	write.table(phe,file=file,quote=F,row.names=F,col.names=T)
	
}