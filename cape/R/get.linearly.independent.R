get.linearly.independent <-
function(geno.matrix){

	# matrixX <- pairscan.obj$geno.for.pairscan
	matrixX <- geno.matrix
	
	if(dim(matrixX)[2] == 1){
		return(matrixX)
		}
	
	#use precision to the 3rd decimal place
	orig.matrixX <- matrixX
	matrixX <- round(matrixX, 3)

	#find the markers without variation
	num.geno <- apply(matrixX, 2, function(x) length(unique(x[!is.na(x)])))
	rejected.markers <- names(num.geno[num.geno < 2])
	good.markers <- num.geno[num.geno >= 2]

	all.pairs <- pair.matrix(names(good.markers))

	if(length(which(is.na(matrixX))) > 0){
		all.cor <- cor(matrixX, use = "pairwise.complete.obs")
		}else{
		all.cor <- cor(matrixX)
		}
	diag(all.cor) <- 0
	
	perfect.cor <- which(abs(all.cor) == 1, arr.ind = TRUE)
	
	#if there are markers with perfect correlation, 
	#remove the first one of the pair
	if(length(perfect.cor) > 0){
		perfect.cor.names <- perfect.cor
		perfect.cor.names[,1] <- colnames(matrixX)[perfect.cor[,1]]
		perfect.cor.names[,2] <- colnames(matrixX)[perfect.cor[,2]]
		rejected.markers <- c(rejected.markers, perfect.cor.names[,1])
		}

	
	rejected.markers <- as.vector(rejected.markers)
	if(length(rejected.markers) > 0){
		rej.markers <- match(sort(unique(rejected.markers)), colnames(matrixX))
		final.geno <- matrixX[,-sort(unique(rej.markers))]
		}else{
			rej.markers <- NULL
			final.geno <- matrixX
			}

	results <- list(final.geno, rej.markers)
	names(results) <- c("independent.markers", "rejected.markers")
	return(results)	

}
