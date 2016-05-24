select.markers.for.pairscan <-
function(data.obj, geno.obj = NULL, singlescan.obj, alpha.thresh = NULL, t.thresh = NULL, specific.markers = NULL, num.markers = NULL, step.size = NULL, tolerance = 10, verbose = TRUE){
	
	if(is.null(alpha.thresh) && is.null(t.thresh) && is.null(specific.markers) && is.null(num.markers)){
		t.thresh = 0
		}

	#take our various parts of the data object
	#for clarity of code later on
	scanone.result <- singlescan.obj$singlescan.results
	
	total.markers <- dim(scanone.result[[1]])[1]
	if(!is.null(num.markers) && num.markers > total.markers){
		num.markers <- total.markers
		}
	
	#we need to use the covariates determined from the 1D scan
	#in the 2D scan, so stop if the 1D scan has not been performed
	if(length(scanone.result) == 0){
		stop("singlescan() must be run before selecting markers for pairscan\n")
		}	

	orig.geno <- get.geno(data.obj, geno.obj)
		

	locus.names <- colnames(orig.geno)
	phenotype.names <- names(scanone.result)
	

	#if we are filtering by the pairs threshold, we 
	#need to find all alleles that exceed the significance
	#threshold for any of the phenotypes. All these alleles
	#get put into a 2D matrix. The rest of the algorithm is
	#identical to cape
	#if we are not thresholding, just flatten the entire
	#array to 2D
			
		
		#============================================================================================
		# Some functions to help parse data structures
		#============================================================================================
		get.sig.allele <- function(result.mat, pairs.thresh){
			sig.alleles <- which(abs(result.mat[,"t.stat"]) >= pairs.thresh, arr.ind = TRUE)
			return(sig.alleles)
			}
		
		
		get.named.alleles <- function(specific.markers){
			marker.locale <- lapply(scanone.result, function(x) which(specific.markers %in% rownames(x)))
			for(i in 1:length(marker.locale)){
				names(marker.locale[[i]]) <- specific.markers
				}
			return(marker.locale)
			}
		
		get.alleles <- function(above.thresh, locus){
			alleles <- c(sort(unique(unlist(lapply(above.thresh, function(x) x[which(x[,"locus"] == locus), "allele"])))))
			allele.list <- as.vector(alleles, mode = "list")
			names(allele.list) <- rep(locus, length(alleles))
			return(allele.list)
			}
					
			
		get.unique.sig <- function(above.thresh){
			all.loci.names <- sort(unique(unlist(sapply(above.thresh, function(x) names(x)))))
			locus.locale <- which(data.obj$marker.num %in% all.loci.names)
			results <- orig.geno[,locus.locale]
			return(results)
			}
								
		#one big function for selecting markers. 
		#this is repeated multiple times if num.markers is specified
		
		get.markers <- function(t.thresh){
			#if no markers have been specified, get above.thresh using the pairs.threshold
			if(is.null(specific.markers)){ 
				if(is.null(t.thresh)){
					pairs.thresh <- 0
					}else{
					pairs.thresh <- t.thresh
					}
					data.obj$t.thresh <- pairs.thresh
					above.thresh <- lapply(scanone.result, function(x) get.sig.allele(x, pairs.thresh))
					if(sum(sapply(above.thresh, length)) == 0 && !is.null(t.thresh)){stop("There are no markers above this threshold. Please choose a lower starting threshold.")}
					
					}else{ 
					#if marker names have been supplied, 
					#create above.thresh using these
					above.thresh <- get.named.alleles(specific.markers)
					}
									
				pair.geno <- get.unique.sig(above.thresh)	
				data.obj$geno.for.pairscan <- pair.geno
			
				#make sure that the genotype matrix is linearly independent
				geno.ind <- get.linearly.independent(pair.geno)
				return(list(geno.ind[[1]], geno.ind[[2]], data.obj))
				}
			
		#============================================================================================
		# end functions			
		#============================================================================================
			
			rejected.markers <- NULL
			
			if(is.null(num.markers)){
				test.obj <- get.markers(t.thresh)
				data.obj <- test.obj[[3]]
				geno <- test.obj[[1]]
				rejected.markers <- test.obj[[2]]
				num.markers.selected <- dim(geno)[2]
				data.obj$geno.for.pairscan <- geno
				if(!is.null(t.thresh)){
					data.obj$t.thresh <- t.thresh
					}
				}else{
				if(verbose){cat("Finding a threshold...\n")}
				if(verbose){cat("thresh\tnum.markers\n")}
				#start with a threshold that should get us near the right number of
				#markers
				scanone.tstats <- unlist(lapply(scanone.result, function(x) x[,3]))
				t.thresh <- as.vector(sort(scanone.tstats, decreasing = TRUE)[num.markers])
				update.val <- 0.01
				less.than <- 1; greater.than <- 1
				iteration <- 1
				test.obj <- get.markers(t.thresh)
				t.thresh <- test.obj[[3]]$t.thresh
				num.markers.selected <- dim(test.obj[[1]])[2]
				if(verbose){cat(t.thresh, "\t", num.markers.selected, "\n")}
				while((num.markers.selected < num.markers - tolerance) || num.markers.selected > num.markers + tolerance){
					if(num.markers.selected < num.markers){
						t.thresh <- t.thresh - update.val
						less.than <- less.than + 1
						}
					if(num.markers.selected > num.markers){
						t.thresh <- t.thresh + update.val
						greater.than <- greater.than + 1
						}
					test.obj <- get.markers(t.thresh)
					num.markers.selected <- dim(test.obj[[1]])[2]
					if(verbose){cat(t.thresh, "  ", num.markers.selected, "\n")}
					iteration <- iteration + 1
					if(greater.than > 2 & less.than > 2){
						stop("I can't get a number of markers within the tolerance. \nPlease adjust the tolerance or the number of markers desired")
						}
					}
				data.obj <- test.obj[[3]]
				geno <- test.obj[[1]]
				rejected.markers <- test.obj[[2]]
				num.markers.selected <- dim(geno)[2]
				data.obj$geno.for.pairscan <- geno
				cat("\n")
				}
			


			if(length(rejected.markers) > 0){
				message("\n", length(rejected.markers), " marker(s) rejected due to linear non-independence.\n For more information see markers.removed.for.non.independence.txt")
				write.table(colnames(geno)[sort(rejected.markers)], "markers.removed.for.non.independence.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
				}
						
			
			
			cat(dim(geno)[2], "markers were selected for the pairscan.\n")
			cat("This makes", choose(dim(geno)[2], 2), "possible pairs.\n")

		return(data.obj)
	
}
