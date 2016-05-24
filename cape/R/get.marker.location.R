get.marker.location <-
function(data.obj, markers){

	is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
	
	if(is.char){
		marker.loc <- data.obj$marker.location[match(markers, data.obj$marker.names)]
		na.locale <- which(is.na(marker.loc))
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			geno.covar <- which(covar.info$covar.type == "g")
			if(length(geno.covar) > 0){
				geno.covar.locale <- match(covar.info$covar.names[geno.covar], data.obj$g.covar[1,])
				geno.covar.loc <- data.obj$g.covar[3,geno.covar.locale]
				geno.covar.pos <- match(data.obj$g.covar[1,], markers)
				marker.loc[geno.covar.pos] <- geno.covar.loc
				}
			}
		na.locale <- which(is.na(marker.loc))
		#if there are still missing values, these are 
		#phenotypic covariates. They get dummy positions
		if(length(na.locale) > 0){ 
			pheno.covar <- which(covar.info$covar.type == "p")
			if(length(pheno.covar) > 0){
				marker.loc[na.locale] <- 1:length(pheno.covar)
				}
			} #end case for if there are still missing markers
		return(marker.loc)
		}
	
	
	if(!is.char){
		marker.loc <- data.obj$marker.location[match(markers, data.obj$marker.num)]
		na.locale <- which(is.na(marker.loc))
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			covar.locale <- which(!is.na(match(colnames(covar.info$covar.table), markers)))
			marker.loc[na.locale] <- covar.info$covar.loc[covar.locale]
			}
		return(marker.loc)
		}
	
}
