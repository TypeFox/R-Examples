filter.maf <-
function(data.obj, geno.obj, maf = 0.05){
		
	geno.mat <- get.geno(data.obj, geno.obj)
	
	removed.markers <- NULL
	
	all.maf <- colMeans(geno.mat, na.rm = TRUE)

	too.high <- which(all.maf > 1-maf)
	if(length(too.high) > 0){
		geno.mat <- geno.mat[,-too.high]
		removed.markers <- c(removed.markers, names(too.high))
		}
		
	too.low <- which(all.maf < maf)
	if(length(too.low) > 0){
		geno.mat <- geno.mat[,-too.low]
		removed.markers <- c(removed.markers, names(too.low))
		}
	
	if(length(removed.markers) > 0){
		cat("Removed", length(removed.markers), "markers due to low minor allele frequency.\n")
		write.table(matrix(sort(removed.markers), ncol = 1), "markers.removed.for.low.MAF.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		#remove the markers with low minor allele frequencey
		data.obj  <- remove.markers(data.obj, removed.markers)
		}
	
	return(data.obj)
	}
