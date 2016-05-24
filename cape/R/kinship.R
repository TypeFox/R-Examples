kinship <-
function(geno, all.chr, chr.which = NULL, chr.pairs = FALSE, sample.kinship = FALSE, num.samples = 100, n.per.sample = 10, run.parallel = TRUE, add.full.kin = FALSE, verbose = FALSE, n.cores = 2){
		
	
	sample.geno <- function(data.step, other.geno){
		# make a sampled similarity matrix for the individuals
		sampled.markers <- sample(1:dim(other.geno)[2], n.per.sample)
		sampled.geno <- other.geno[,sampled.markers]
		return(sampled.geno)
		}
	
	
	sample.K <- function(geno.sample, use.cov = FALSE){
		if(!use.cov){
			#This method seems to be faster with small matrices
			K.sample <- (geno.sample %*% t(geno.sample))/ncol(geno.sample)
			}else{
			K.sample <- cov(t(geno.sample)) + (colMeans(t(geno.sample)) %*% t(colMeans(t(geno.sample))))
			}
		return(K.sample)
		}


	combine.samples <- function(sample.list){
		if(run.parallel){
			cl <- makeCluster(n.cores)
			registerDoParallel(cl)
			total.K <- foreach(n = sample.list, .combine = "+") %dopar% {
				sample.K(n)
				}
			stopCluster(cl)
			}else{
			total.K <- matrix(0, nrow = nrow(geno), ncol = nrow(geno))
			for(n in 1:num.samples){
				if(verbose){report.progress(n, num.samples)}
				total.K <- total.K + sample.K(sample.list[[n]])
				}
			}		
			final.K <- total.K/num.samples
			return(final.K)	
			}
	

	get.g=function(pair = NULL){
		#set pair to null if we are only looking at the covariate
		if(!is.null(pair)){
			#find the locations of the markers not in the pair
			not.pair.chr <- which(!all.chr %in% pair)
			other.geno = geno[,not.pair.chr,drop=FALSE]
			}else{
			other.geno <- geno	
			}
		
		if(sample.kinship){
			if(verbose){cat("\nSampling genotypes...\n")}
			geno.samples <- lapply(1:num.samples, function(x) sample.geno(x, other.geno))
			K <- combine.samples(geno.samples)
			}else{
			K <- cov(t(other.geno)) + (colMeans(t(other.geno)) %*% t(colMeans(t(other.geno))))				
			}
		
		return(K)
		}
	
	if(is.null(chr.which)){
		chr <- unique(all.chr)
		}else{
		chr <- chr.which	
		}
	
	if(chr.pairs){
		pairs.mat <- pair.matrix(chr, self.pairs = TRUE)
		}else{
		pairs.mat <- cbind(chr, chr)	
		}

	list.names <- apply(pairs.mat, 1, function(x) paste(x, collapse = ","))
	
	kin.list <- vector(mode = "list", length = dim(pairs.mat)[1])
	names(kin.list) <- list.names

	if(verbose){cat("Calculating kinship matrices...\n")}
	for(k in 1:dim(pairs.mat)[1]){
		if(verbose){cat("\nChromosome(s)", unique(pairs.mat[k,], "\n"))}
		kin.list[[k]] <- get.g(pair = pairs.mat[k,])
		}
		if(verbose){cat("\n")}

	#if there are covariates, also calulate the kinship matrix for the full genome
	if(add.full.kin){
		if(verbose){cat("\nCalculating kinship matrix with full genome for testing covariates.\n")}
		kin.list[[k+1]] <- get.g(pair = NULL)
		names(kin.list)[k+1] <- "full.geno"
		}
	
	
	return(kin.list)
	
}
