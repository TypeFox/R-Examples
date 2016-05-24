get.geno <-
function(data.obj, geno.obj){
		
	
	if(missing(geno.obj) || is.null(geno.obj)){
		geno <- data.obj$geno
		geno.names <- dimnames(geno)
		}else{
		geno <- geno.obj$geno
		geno.names <- list(rownames(data.obj$pheno), data.obj$marker.names)
		}
		
	
	if(is.null(geno)){
		stop("I can't find the genotype data. Please make sure it is in either data.obj or geno.obj.")
		}
	
	ind.dim <- 1
	locus.dim <- 2

	#subset the genotype object to match the 
	#individuals and markers we want to scan

	ind.locale <- match(geno.names[[ind.dim]], dimnames(geno)[[ind.dim]])
	locus.locale <- match(geno.names[[locus.dim]], dimnames(geno)[[locus.dim]])
		
	gene <- geno[ind.locale, locus.locale]
			
	return(gene)	
}
