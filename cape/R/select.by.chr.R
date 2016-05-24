select.by.chr <-
function(data.obj, chr){
	
	geno <- data.obj$geno
	
	chr.list <- data.obj$chromosome
	chr.locale <- which(chr.list %in% chr)
	
	sub.geno <- geno[,chr.locale]
	data.obj$geno <- sub.geno
	data.obj$chromosome <- data.obj$chromosome[chr.locale]
	data.obj$marker.location <- data.obj$marker.location[chr.locale]
	data.obj$marker.names <- data.obj$marker.names[chr.locale]

	return(data.obj)
	
	}
