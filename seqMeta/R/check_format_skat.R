check_format_skat <- function(Z, SNPInfo, mod0, aggregateBy, snpNames, formula) {
	if(length(stats::residuals(mod0)) != nrow(Z)) stop("Number of genotypes is not equal to number of phenotypes!")
	
	if(!is.null(stats::na.action(mod0))){ 
	  stop(paste0("Some observations in '", 
	              utils::capture.output(print(formula)), 
	              "' are missing...\n Complete data in the null model is required. Please remove, and subset genotypes accordingly"))
	}
  
	percent_miss <- mean(colnames(Z) %in% SNPInfo[,snpNames])
	if(percent_miss == 0) stop("Column names of Z must correspond to 'snpNames' field in SNPInfo file!")
	if(percent_miss < 1) warning(paste((1-percent_miss)*ncol(Z), "snps are not in SNPInfo file!"))
	
	if( min(Z,na.rm=TRUE) < 0 || max(Z,na.rm=TRUE) > 2) {
    warning("Genotype possibly not dosage matrix! Alleles should be coded 0/1/2")
	}
	if(anyNA(Z)) {
    warning("Some missing genotypes - will be imputed to average dose")
	}
}
