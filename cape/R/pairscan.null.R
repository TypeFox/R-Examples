pairscan.null <-
function(data.obj, geno.mat, covar = NULL, scan.what = c("eigentraits", "raw.traits"), total.perm = NULL, n.top.markers = 50, max.pair.cor, min.per.genotype, verbose = FALSE, n.cores = NULL){

	covar.table = NULL #for appeasing R CMD check


	pairscan.obj <- vector(mode = "list", length = 2)
	names(pairscan.obj) <- c("pairscan.perm", "pairs.tested.perm")

	if(is.null(total.perm)){
		stop("The total number of permutations must be specified.")
		}

	pheno <- get.pheno(data.obj, scan.what)
	
	pheno.names <- colnames(pheno)
	num.pheno <- dim(pheno)[2]
	#use the full genotype matrix to select 
	#markers for generating the null in the 
	#pairscan
	
	
	covar.info <- get.covar(data.obj, covar)
	for(i in 1:length(covar.info)){
		assign(names(covar.info)[i], covar.info[[i]])
		}

	#make a list to hold the results. 
	#Tables from each of the phenotypes will be
	#put into this list
	results.perm.list <- vector(mode = "list", length = num.pheno)
 	names(results.perm.list) <- colnames(pheno)
 	
	#generate the null distribution for the pairscan 
	#do a singlescan on each permuted trait.
	#find the top n markers
	#combine these into a unique list
	#do the pairscan on these markers for all traits

	if(verbose){
		cat("\nGenerating null distribution...\n")
		}

	all.pairs.tested <- NULL
	
	null.markers <- matrix(NA, ncol = 1, nrow = n.top.markers*num.pheno)
	final.perm <- 1
	while(final.perm < total.perm){
		start.marker = 1
		perm.order <- sample(1:dim(pheno)[1])
		for(p in 1:num.pheno){ 
			if(n.top.markers < dim(geno.mat)[2]){
				if(verbose){cat("\tSingle-marker scan of permuted", colnames(pheno)[p], "...\n")}
				single.scan.result <- one.singlescan(phenotype.vector = pheno[perm.order,p], genotype.mat = geno.mat, covar.vector = covar.table, n.cores = n.cores)
			
				#find the loci with the largest effect				
				geno.order <- order(single.scan.result[,"t.stat"], decreasing = TRUE)
				#record the loci found for each phenotype
				null.markers[start.marker:(start.marker+n.top.markers-1),1] <- rownames(single.scan.result)[geno.order[1:n.top.markers]]

				}else{
				null.markers[start.marker:(start.marker+n.top.markers-1),1] <- colnames(geno.mat)
				}
			
					

			start.marker = start.marker + n.top.markers
			}

				#find the unique top markers
				u_markers <- unique(null.markers[(which(!is.na(null.markers)))])
			
				
				#make a matrix with the genotype of each max locus in the top n loci
				#this will be a 2D matrix with a row for each individual, and a column 
				#for each max locus
				top.geno.mat <- matrix(unlist(lapply(1:length(u_markers), function(x) geno.mat[,u_markers[x]])), ncol = length(u_markers), byrow = FALSE)				
				colnames(top.geno.mat) <- u_markers
				rownames(top.geno.mat) <- rownames(geno.mat)
		
				if(verbose){cat("\tGetting markers for permuted pairscan...\n")}
				top.marker.pairs <- get.pairs.for.pairscan(geno = top.geno.mat, min.per.genotype = min.per.genotype, max.pair.cor = max.pair.cor, verbose = verbose, n.cores = n.cores)

				#we don't need to do a million extra permutations
				#so trim the final pair matrix down to get only
				#the specified number of permutations
				if(final.perm+dim(top.marker.pairs)[1] > total.perm){
					num.needed <- total.perm - final.perm
					top.marker.pairs <- top.marker.pairs[sample(1:dim(top.marker.pairs)[1], num.needed),,drop=FALSE]
					}

				if(verbose){cat("\tTesting", dim(top.marker.pairs)[1], "pairs...\n")}
				all.pairs.tested <- rbind(all.pairs.tested, top.marker.pairs)
				
			#run the pairscan for each permuted phenotype and the pairs we just found
			for(p in 1:num.pheno){
				if(verbose){cat("\tMarker-pair scan of permuted", colnames(pheno)[p], "...\n")}
				#run a pairscan on these markers and each permuted phenotype
				pairscan.results <- one.pairscan.parallel(data.obj, phenotype.vector = pheno[perm.order,p], genotype.matrix = top.geno.mat, paired.markers = top.marker.pairs, n.perm = 0, verbose = FALSE, covar.vector = NULL, n.cores = n.cores)
				
				#integrate the results into the permutation object
				one.perm <- pairscan.results[[1]]
				if(final.perm == 1){ #if this is the first time through, just copy the results into the results.perm.list
					results.perm.list[[p]] <- one.perm
					}else{
					for(i in 1:length(one.perm)){
						results.perm.list[[p]][[i]] <- rbind(results.perm.list[[p]][[i]], one.perm[[i]])
						}
					}
				}
				final.perm <- dim(results.perm.list[[1]][[1]])[1] #end looping through phenotypes
				if(verbose){cat("\t", final.perm, " permutations: ", round((final.perm/total.perm)*100), "%\n", sep = "")} 
				} #end when we have enough permutations
		
		
	names(results.perm.list) <- pheno.names
	pairscan.obj$pairscan.perm <- results.perm.list
	pairscan.obj$pairs.tested.perm <- all.pairs.tested
	return(pairscan.obj) #and return it
	
}
