get.marker.name <-
function(data.obj, markers){
	
	is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))

	if(is.char){
		#use the marker.names vector first to translate
		marker.name <- data.obj$marker.names[match(markers, data.obj$marker.names)]
		
		#if there are any markers we didn't translate, look in the 
		#covariate tables for marker numbers
		na.locale <- which(is.na(marker.name))
		
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			marker.name[na.locale] <- covar.info$covar.names[match(markers[na.locale], covar.info$covar.names)]
			}
	}else{
		#use the marker.names vector first to translate
		marker.name <- data.obj$marker.names[match(markers, data.obj$marker.num)]
		
		#if there are any markers we didn't translate, look in the 
		#covariate tables for marker numbers
		na.locale <- which(is.na(marker.name))
		
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			marker.name[na.locale] <- covar.info$covar.names[match(markers[na.locale], colnames(covar.info$covar.table))]
			}
		}

	return(marker.name)

	
	
}
