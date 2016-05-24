get.marker.chr <-
function(data.obj, markers){
	
	is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
	
	if(is.char){
		marker.chr <- data.obj$chromosome[match(markers, data.obj$marker.names)]
		na.locale <- which(is.na(marker.chr))
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			geno.covar <- which(covar.info$covar.type == "g")
			if(length(geno.covar) > 0){
				unassigned <- markers[na.locale]
				unassigned.locale <- match(unassigned, data.obj$g.covar[1,])
				geno.covar.chr <- data.obj$g.covar[2,unassigned.locale]
				marker.chr[na.locale] <- geno.covar.chr
				}
			}
		#The rest are phenotypic covariates which are assigned chr 0
		na.locale <- which(is.na(marker.chr))
		if(length(na.locale) > 0){
			marker.chr[na.locale] <- 0
			}
		return(marker.chr)
		}
	
	
	
	if(!is.char){
		marker.chr <- data.obj$chromosome[match(markers, data.obj$marker.num)]
		na.locale <- which(is.na(marker.chr))
		if(length(na.locale) > 0){
			covar.info <- get.covar(data.obj)
			geno.covar <- which(covar.info$covar.type == "g")
			if(length(geno.covar) > 0){
				unassigned <- markers[na.locale]
				unassigned.locale <- match(unassigned, colnames(data.obj$g.covar))
				geno.covar.chr <- data.obj$g.covar[2,unassigned.locale]
				marker.chr[na.locale] <- geno.covar.chr
				}
			}
		#The rest are phenotypic covariates which are assigned chr 0
		na.locale <- which(is.na(marker.chr))
		if(length(na.locale) > 0){
			marker.chr[na.locale] <- 0
			}
		return(marker.chr)
		}






}
