writePopulation <-
function(data.obj, filename, geno.to.convert = NULL, conversion = NULL){
	

	geno <- data.obj$geno
	pheno <- data.obj$pheno
	chr <- data.obj$chromosome
	loc <- data.obj$marker.location
	marker.names <- data.obj$marker.names

	if(!is.null(geno.to.convert)){
		if(length(geno.to.convert) != length(conversion)){
			stop("The geno.to.convert and conversion vectors must be the same length.")
			}
			for(i in 1:length(geno.to.convert)){
				geno.locale <- which(geno == geno.to.convert[i])
				geno[geno.locale] <- conversion[i]
				}	
			}

	#build a matrix out of the object so we can write it out to a csv file
	final.geno <- rbind(chr, loc, geno)
	final.geno[which(is.na(final.geno))] <- "-"
	colnames(final.geno) <- marker.names
	
	pheno[which(is.na(pheno))] <- "-"
	pheno.padding <- matrix(NA, nrow = 2, ncol = dim(pheno)[2])
	final.pheno <- rbind(pheno.padding, pheno)
	
	
	final.data <- cbind(final.pheno, final.geno)
	write.table(final.data, file = filename, sep = ",", quote = FALSE, row.names = FALSE, na = "")

	
	
}
