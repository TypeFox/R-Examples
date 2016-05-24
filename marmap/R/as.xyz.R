as.xyz <- function(bathy) {

	if (class(bathy)!="bathy") stop("Objet is not of class bathy")
	
	lon <- as.numeric(rownames(bathy))
	lat <- as.numeric(colnames(bathy))
	
	xyz <- data.frame(expand.grid(lon,lat),as.vector(bathy))
	xyz <- xyz[order(xyz[,2], decreasing=TRUE),]
	names(xyz) <- c("V1","V2","V3")
	rownames(xyz) <- 1:nrow(xyz)

	return(xyz)
}