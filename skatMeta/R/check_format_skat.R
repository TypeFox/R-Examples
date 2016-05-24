check_format_skat <-
function(Z, SNPInfo, mod0, aggregateBy, snpNames){
	if(length(residuals(mod0)) != nrow(Z)) stop("Number of genotypes is not equal to number of phenotypes!")
	
	percent_miss <- mean(colnames(Z) %in% SNPInfo[,snpNames])
	if(percent_miss == 0) stop("Column names of Z must correspond to 'snpNames' field in SNPInfo file!")
	if(percent_miss < 1) warning(paste((1-percent_miss)*ncol(Z), "snps are not in SNPInfo file!"))
	
	if({rr <- range(Z,na.rm=TRUE); rr[1] < 0 & rr[2] > 2}) warning("Genotype possibly not dosage matrix! Alleles should be coded 0/1/2")
	if(any(is.na(Z))) warning("Some missing genotypes - will be imputed to average dose")
}
