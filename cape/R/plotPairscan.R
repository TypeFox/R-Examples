plotPairscan <-
function(data.obj, pairscan.obj, phenotype = NULL, standardized = TRUE, show.marker.labels = FALSE, show.chr = TRUE, label.chr = TRUE, pdf.label = "Pairscan.Regression.pdf", neg.col = "blue", pos.col = "brown", col.pal = "dark", verbose = FALSE){
	
		
	#get the markers used in the pair scan and sort them.
	markers.used <- sort(as.numeric(unique(c(pairscan.obj$pairscan.results[[1]][[1]][,1], pairscan.obj$pairscan.results[[1]][[1]][,2]))))
	
	#get coordinates of the chromosome boundaries
	if(show.chr){
		chromosomes <- get.marker.chr(data.obj, markers.used)
		u_chr <- unique(chromosomes)
		chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
		chr.boundaries <- c(0, chr.boundaries)
		if(label.chr){
			chr.names <- unique(chromosomes)
			chr.names[which(chr.names == 0)] <- "c"
			}else{
			chr.names <- NULL	
			}
		}else{
		chr.boundaries <- NULL
		chr.names <- NULL
		}

	pairscan.result <- pairscan.obj$pairscan.results
	
	if(is.null(pairscan.result)){
		stop("pairscan() must be run before plotPairscan()")
		}
	
	if(is.null(phenotype)){
		phenotype <- names(pairscan.result)
		}
		
	pheno.num <- which(names(pairscan.result) %in% phenotype)
	
	if(length(pheno.num) < length(phenotype)){
		not.found <- which(!(phenotype %in% names(pairscan.result)))
		message("I couldn't find the following phenotypes:")
		cat(phenotype[not.found], sep = "\n")
		stop()
		}
		
	num.pheno <- length(pheno.num)
	
	marker.pairs <- pairscan.obj$pairscan.results[[1]][[1]][,1:2]

		
	#collect the results, so we can put them on the same scale
	all.results.mats <- list()
	min.x <- 0
	max.x <- 0
	#for each phenotype scanned
	for(p in pheno.num){
		if(verbose){cat("\n", phenotype[p], ":\n", sep = "")}
		#build a results matrix
		results.mat <- matrix(0, length(markers.used), length(markers.used))
		colnames(results.mat) <- rownames(results.mat) <- markers.used
		
		if(standardized){
			pair.int <- as.numeric(pairscan.result[[p]][[1]][, "marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][,"marker1:marker2"])
			}else{
			pair.int <- as.numeric(pairscan.result[[p]][[1]][,"marker1:marker2"])
			}
		
		for(i in 1:length(markers.used)){
			m1.locale <- which(pairscan.result[[p]][[1]][,1] == markers.used[i])
			m2.locale <- which(pairscan.result[[p]][[1]][,2] == markers.used[i])
			if(length(m1.locale) > 0){
				results.mat[unique(as.character(pairscan.result[[p]][[1]][m1.locale,1])), as.character(pairscan.result[[p]][[1]][m1.locale,2])] <- pair.int[m1.locale]
				}
			if(length(m2.locale) > 0){
				results.mat[as.character(pairscan.result[[p]][[1]][m2.locale,1]), unique(as.character(pairscan.result[[p]][[1]][m2.locale,2]))] <- pair.int[m2.locale]
				}
			}
	
	
		diag(results.mat) <- NA
		all.results.mats[[p]] <- results.mat
		min.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)*-1
		max.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)
		}	
	


		pdf(pdf.label, width = 7, height = 6)
		for(p in 1:length(pheno.num)){
			results <- all.results.mats[[pheno.num[p]]] + t(all.results.mats[[pheno.num[p]]])
			myImagePlot(x = results, xlab = "marker1", ylab = "marker2", main = phenotype[p], min.x = min.x, max.x = max.x, show.labels = show.marker.labels, chromosome.coordinates = chr.boundaries, chr.names = chr.names, neg.col = neg.col, pos.col = pos.col, col.pal = col.pal)
			}
		dev.off()
		
		
	
	
}
