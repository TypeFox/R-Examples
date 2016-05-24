get.geno.with.covar <-
function(data.obj, geno.obj = NULL, g.covar = TRUE, p.covar = TRUE, for.pairscan = TRUE){
	
	
		covar.info <- get.covar(data.obj, covar = NULL)
				
		covar.locale <- NULL
		if(g.covar){
			covar.locale <- c(covar.locale, which(covar.info$covar.type == "g"))
			}
		if(p.covar){
			covar.locale <- c(covar.locale, which(covar.info$covar.type == "p"))
			}

		
		if(for.pairscan){
			geno <- data.obj$geno.for.pairscan			
			}else{
			geno <- get.geno(data.obj, geno.obj)
			}
			
		is.char <- as.logical(is.na(suppressWarnings(as.numeric(colnames(geno)[1]))))
		if(is.char){
			colnames(covar.info$covar.table) <- covar.info$covar.names
			}

		geno <- cbind(geno, covar.info$covar.table[,covar.locale,drop=FALSE])

		#if there are marker covariates make sure these
		#are placed in the right order in the genotype
		#matrix
		if(g.covar && length(which(covar.info$covar.type == "g"))  > 0){
			new.geno.chr <- get.marker.chr(data.obj, colnames(geno))
			new.geno.pos <- get.marker.location(data.obj, colnames(geno))
			marker.pos.table <- cbind(new.geno.chr, new.geno.pos)
			marker.pos.table <- sortByThenBy(marker.pos.table, col.type = c("n", "n"), return.order = TRUE)
			for(i in 1:2){
				geno <- geno[,marker.pos.table[,i]]
				}

			#now we need to make sure any phenotypic covariates
			#are at the end
			new.geno.chr <- get.marker.chr(data.obj, colnames(geno))
			pheno.covar.locale <- which(new.geno.chr == 0)
			if(length(pheno.covar.locale) > 0){
				geno <- cbind(geno[,-pheno.covar.locale,drop=FALSE],geno[,pheno.covar.locale,drop=FALSE])
				}
			}	

		return(geno)
	
}
