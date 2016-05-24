sortCross <-
function(data.obj, geno.obj){
	
	geno <- get.geno(data.obj, geno.obj)
	new.geno <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
	new.chr <- rep(NA, length(data.obj$chromosome))
	new.marker.location <- rep(NA, length(data.obj$marker.location))

	
	u_chr <- sort(unique(data.obj$chromosome))
	if(u_chr[1] == 0){
		u_chr <- c(u_chr[-1], 0)
		}
		
	for(i in 1:length(u_chr)){
		chr.locale <- which(data.obj$chromosome == u_chr[i])
		marker.order <- order(data.obj$marker.location[chr.locale])
		new.chr[chr.locale] <- data.obj$chromosome[chr.locale[marker.order]]
		new.marker.location[chr.locale] <- data.obj$marker.location[chr.locale[marker.order]]
		new.geno[,chr.locale] <- geno[,chr.locale[marker.order]]
		}
	
	colnames(new.geno) <- 1:dim(new.geno)[2]
	rownames(new.geno) <- 1:dim(new.geno)[1]
	data.obj$geno <- new.geno
	data.obj$marker.location <- new.marker.location
	return(data.obj)
	
}
