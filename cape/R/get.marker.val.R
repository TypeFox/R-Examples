get.marker.val <-
function(data.obj, geno.obj = NULL, markers){
	
	geno <- get.geno(data.obj, geno.obj)
	
	is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))	

	marker.vals <- matrix(NA, nrow = nrow(geno), ncol = length(markers))
	
	if(is.char){
		#use the markers vector first to translate
		marker.num <- data.obj$marker.num[match(markers, data.obj$marker.names)]
		
		#grab values for the markers we found
		not.na.locale <- which(!is.na(marker.num))
		marker.idx <- get.marker.idx(data.obj, marker.num[not.na.locale])
		marker.vals[,not.na.locale] <- geno[,marker.idx]

		#if there are any markers we didn't translate, look in the 
		#covariate tables for marker numbers
		na.locale <- which(is.na(marker.num))
		
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			marker.vals[,na.locale] <- covar.info$covar.table[,match(markers[na.locale], covar.info$covar.names)]
			}
		}else{
		#use the markers vector first to translate
		marker.num <- data.obj$marker.num[match(markers, data.obj$marker.num)]
		
		#grab values for the markers we found
		not.na.locale <- which(!is.na(marker.num))
		marker.idx <- get.marker.idx(data.obj, marker.num[not.na.locale])
		marker.vals[,not.na.locale] <- geno[,marker.idx]
		
		#if there are any markers we didn't translate, look in the 
		#covariate tables for marker numbers
		na.locale <- which(is.na(marker.num))
		
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			marker.vals[,na.locale] <- covar.info$covar.table[,match(markers[na.locale], colnames(covar.info$covar.table))]
			}
		}

		colnames(marker.vals) <- markers
		return(marker.vals)
	
	
}
