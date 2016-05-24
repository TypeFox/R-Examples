get.marker.num <-
function(data.obj, markers){

	is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
	
	
	if(is.char){
		#use the markers vector first to translate
		marker.num <- data.obj$marker.num[match(markers, data.obj$marker.names)]
		
		#if there are any markers we didn't translate, look in the 
		#covariate tables for marker numbers
		na.locale <- which(is.na(marker.num))
		
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			marker.num[na.locale] <- colnames(covar.info$covar.table)[match(markers[na.locale], covar.info$covar.names)]
			}
		}else{
		#use the markers vector first to translate
		marker.num <- data.obj$marker.num[match(markers, data.obj$marker.num)]
		
		#if there are any markers we didn't translate, look in the 
		#covariate tables for marker numbers
		na.locale <- which(is.na(marker.num))
		
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			marker.num[na.locale] <- colnames(covar.info$covar.table)[match(markers[na.locale], colnames(covar.info$covar.table))]
			}
		}


	return(marker.num)
}
