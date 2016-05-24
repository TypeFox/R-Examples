genome.wide.threshold.1D.parallel <-
function(data.obj, geno.mat, n.perm = 1000, alpha = c(0.01, 0.05), scan.what = c("eigentraits", "raw.traits"), verbose = FALSE, n.cores = 2){
	
	#for testing
	# gene <- geno.mat[,1:100]

	gene <- geno.mat
	
	#calculate the numbers of markers, phenotypes and samples
	n.gene <- dim(gene)[2]


	pheno <- get.pheno(data.obj, scan.what)	
	num.samples <- dim(pheno)[1]
	n.phe <- dim(pheno)[2]
	
	#====================================================
	# internal functions
	#====================================================
	get.stat <- function(regression){
			if(dim(summary(regression)$coefficients)[1] == 2){
				stat <- summary(regression)$coefficients[2,1]/summary(regression)$coefficients[2,2]
				}else{
					stat <- NA
					}
				return(stat)
			}
			
	get.s <- function(evd.result, alpha){
		s <- qgev(1-alpha,loc=evd.result$estimate[1], scale=evd.result$estimate[2], shape=evd.result$estimate[3], lower.tail = TRUE)
		return(s)
		}

	#====================================================
		
	 	cl <- makeCluster(n.cores)
	 	registerDoParallel(cl)
	 	
		perm.max <- matrix(NA, ncol = n.phe, nrow = n.perm)
		for(n in 1:n.perm){
			if(verbose){report.progress(n, n.perm)}
			sampled.vector <- sample(1:num.samples)
			pheno.perm <- pheno[sampled.vector,]
			pheno.perm.mat <- sampled.vector

	 		for(et in 1:dim(pheno)[2]){
	  			regress.list <- foreach(g = gene, .combine = "c") %dopar% {
	  				get.stat(lm(pheno.perm[,et]~g))
	  				}
			  	perm.max[n,et] <- max(regress.list, na.rm = TRUE)
  				}
			}      
		stopCluster(cl)

        
	#apply the extreme value distribution to the results
	evd <- apply(perm.max, 2, function(x) fgev(x, std.err = FALSE))

	
	s <- vector(mode = "list", length = length(alpha))
	for(a in 1:length(alpha)){
		s[[a]] <- as.vector(sapply(evd, function(x) get.s(x, alpha[a])))
		}
	
	#calculate one threshold over all phenotypes
	thresholds <- lapply(s, mean)
	names(thresholds) <- alpha
		
	if(verbose){
		cat("\n") #make sure the prompt is on the next line at the end of everything
		}
		
	return(thresholds)


}
