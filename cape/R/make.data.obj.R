make.data.obj <-
function(pheno.obj, geno.obj){
	
	data.obj <- list("pheno" = pheno.obj)
	
	data.obj$marker.names <- geno.obj$marker.names
	data.obj$marker.num <- 1:length(geno.obj$marker.names)
	data.obj$chromosome <- geno.obj$chromosome
	data.obj$marker.location <- geno.obj$marker.location

	pheno.ind <- rownames(data.obj$pheno)
	geno.ind <- rownames(geno.obj$geno)
	
	ind.locale <- match(geno.ind, pheno.ind)
	# cbind(pheno.ind[ind.locale], geno.ind)
	ordered.pheno <- pheno.obj[ind.locale[which(!is.na(ind.locale))],]

	data.obj$pheno <- ordered.pheno
	
	return(data.obj)

	
}
