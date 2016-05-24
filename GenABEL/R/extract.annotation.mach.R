#' extracts SNP annotation from MACH/HapMap legend files
#' 
#' This function extracts SNP annotation from MACH info and HapMap legend files. 
#' 
#' @param infofile MACH (ml)info-file name
#' @param legendfile HapMap legend file name
#' @param chromosome chromosome
#' 
#' @return data frame containing annotaton
#' 
#' @author Yurii Aulchenko
#' 
#' @keywords IO manip
#' 

extract.annotation.mach <- function(infofile,legendfile,chromosome=NA)
{
	
	if (!file.exists(legendfile)) stop(paste("legendfile",legendfile,"does not exist"))
	if (!file.exists(infofile)) stop(paste("infofile",infofile,"does not exist"))
	
	# SNP	Al1	Al2	Freq1	MAF	Quality	Rsq
	tmp <- extract_text_file_columns(infofile,c(1:7))[-1,]
	infoannot <- as.data.frame(tmp,stringsAsFactors=FALSE)
	names(infoannot) <- c("name","A1","A0","Freq1","MAF","Quality","Rsq")
	
	
	#rs	position	0	1
	tmp <- extract_text_file_columns(legendfile,c(1:2))[-1,]
	legendannot <- as.data.frame(tmp,stringsAsFactors=FALSE)
	rm(tmp);gc()
	names(legendannot) <- c("name","pos")
	# merge annotation
	
	ord <- order(as.character(infoannot$name))
	
	mrg <- merge(infoannot,legendannot,by="name",all.x=TRUE,all.y=F)
	rownames(mrg) <- mrg$name
	mrg <- mrg[as.character(infoannot$name),]
	class(mrg$pos) <- class(mrg$Freq1) <- class(mrg$Rsq) <- class(mrg$Quality) <- "numeric"
	if (!missing(chromosome)) mrg$chr <- chromosome
	return(mrg)
	
}