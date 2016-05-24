remove.markers <-
function(data.obj, markers.which){
		
			
	marker.locale <- get.marker.idx(data.obj, markers.which)


	geno.locale <- which(names(data.obj) == "geno")
	if(length(geno.locale) > 0){
		data.obj[[geno.locale]] <- data.obj[[geno.locale]][,-marker.locale,drop=FALSE]
		}
	data.obj$marker.num <- data.obj$marker.num[-marker.locale]
	data.obj$marker.names <- data.obj$marker.names[-marker.locale]
	data.obj$chromosome <- data.obj$chromosome[-marker.locale]
	data.obj$marker.location <- data.obj$marker.location[-marker.locale]



	if(!is.null(data.obj$geno.for.pairscan)){
		is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers.which[1]))))
		if(!is.char){
			marker.names <- get.marker.name(data.obj, markers.which)
			pair.geno.locale <- which(colnames(data.obj$geno.for.pairscan) %in% marker.names)
			}else{
			marker.num <- get.marker.num(data.obj, markers.which)
			pair.geno.locale <- which(colnames(data.obj$geno.for.pairscan) %in% markers.which)
			}

	#remove markers from the geno.for.pairscan table if they are there
		if(length(pair.geno.locale) > 0){
			data.obj$geno.for.pairscan <- data.obj$geno.for.pairscan[,-pair.geno.locale]
			}
		}
		
	return(data.obj)
	}
