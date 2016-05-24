#' extracts SNP annotation from IMPUTE files
#' 
#' This function extracts SNP annotation information from IMPUTE files. 
#' The major problem at the moment that info-file format of IMPUTE is 
#' a little bit unstable (reported information and column order 
#' varies between impute v1, v2, and beta-version). Therefore take 
#' special care to read specification of 'order_info_snp_pos_freq1_info_qual_type'
#' 
#' 
#' @param genofile IMPUTE genotype file name
#' @param infofile IMPUTE info-file name
#' @param chromosome chromosome
#' @param order_geno_snp_a0_a1 which columns to extract from geno-file, 
#' and what is the order for snp name, a0, and a0? (default is OK)
#' @param skip_geno how many lines of geno-file are to 
#' be skept? (default is OK)
#' @param order_info_snp_pos_freq1_info_qual_type 
#' which columns to extract from info-file, 
#' and what is the order for SNP name, position, frequency of allele 1, 
#' info (Rsq), and quality (average max post prob)? Dafult works for 
#' IMPUTE v2.0, but has to be changed for other versions. Always check!
#' @param skip_info how many lines of info-file are to be skept 
#' before information starts? IMPUTE v2.0 has a header line, therefore 
#' skip_info=1 works fine; this may be different for other versions 
#' of IMPUTE
#' @param allow_duplicated_names if duplicated SNP names are allowed
#' (same order in geno and info- files is assumed then)
#' 
#' @return data frame containing annotaton
#' 
#' @author Yurii Aulchenko
#' 
#' @keywords IO manip
#' 

extract.annotation.impute <- function(genofile,infofile,chromosome=NA,
		order_geno_snp_a0_a1=c(2,4:5),skip_geno=0,
		order_info_snp_pos_freq1_info_qual_type=c(2:7),skip_info=1,
		allow_duplicated_names=FALSE
)
{
	if (!require(DatABEL))
		stop("this function requires DatABEL package to be installed")
	
	if (!file.exists(genofile)) stop(paste("genofile",genofile,"does not exist"))
	if (!file.exists(infofile)) warning(paste("infofile",infofile,"does not exist; skipping that info"))
	
	# rsname position Al1 Al2
	if (skip_geno>0)
		tmp <- extract_text_file_columns(genofile,order_geno_snp_a0_a1)[-c(1:skip_geno),]
	else
		tmp <- extract_text_file_columns(genofile,order_geno_snp_a0_a1)
	genoannot <- as.data.frame(tmp,stringsAsFactors=FALSE)
	names(genoannot) <- c("name","A0","A1")
	if (!file.exists(infofile)) {
		return(genoannot)
	}
	#1snp_id 2rs_id 3position 4exp_freq_a1 5info 6certainty 7type info_type1 concord_type1 r2_type1 info_type0 concord_type0 r2_type0
	#rs_id position exp_freq_a1 info certainty type
	if (skip_info>0)
		tmp <- extract_text_file_columns(infofile,order_info_snp_pos_freq1_info_qual_type)[-c(1:skip_info),]
	else 
		tmp <- extract_text_file_columns(infofile,order_info_snp_pos_freq1_info_qual_type)
	
	infoannot <- as.data.frame(tmp,stringsAsFactors=FALSE)
	rm(tmp);gc()
	names(infoannot) <- c("name","pos","Freq1","Rsq","Quality","type")
	#print(infoannot[1:5,])
	# merge annotation
	if (!all(genoannot$name %in% infoannot$name)) {
		warning("not all snps are in info-file")
		warning(paste("missing snps:",genoannot$name[which(!(genoannot$name %in% infoannot$name))]))
	}
	if (allow_duplicated_names) {
		if (dim(genoannot)[1] != dim(infoannot)[1]) 
			stop("number of SNPs in geno- and info- files different")
		mrg <- cbind(genoannot,infoannot[,c(2:dim(infoannot)[2])])
	} else {
		mrg <- merge(genoannot,infoannot,by="name",all.x=TRUE,all.y=F)
		rownames(mrg) <- mrg$name
		mrg <- mrg[as.character(genoannot$name),]
	}
	class(mrg$pos) <- class(mrg$Freq1) <- class(mrg$Rsq) <- class(mrg$Quality) <- "numeric"
	if (!missing(chromosome)) mrg$chr <- chromosome
	return(mrg)
}