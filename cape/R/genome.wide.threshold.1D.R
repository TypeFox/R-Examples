genome.wide.threshold.1D <-
function(data.obj, geno.mat, n.perm = 1000, alpha = c(0.01, 0.05), scan.what = c("eigentraits", "raw.traits"), verbose = FALSE){
	
	# require("evd")
	gene <- geno.mat
	#calculate the numbers of markers, phenotypes and samples
	n.gene <- dim(gene)[2]


	#pull out genotype and phenotype data for
	#clarity of code later.
	#If the user has not specified a scan.what,
	#from the outer function (singlescan.R),
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentrais, otherwise, use raw phenotypes
		type.choice <- c(grep("eigen", scan.what), grep("ET", scan.what), grep("et", scan.what))
		if(length(type.choice) > 0){
			pheno <- data.obj$ET
			}else{
				pheno <- data.obj$pheno
				}
	
		num.samples <- dim(pheno)[1]
	
	
	#make additional objects to be used in generating a null for the pairscan
	#these include a matrix containing each permutation of the phenotype
	#and a matrix of the rank of each marker in terms of standardized effect size
	pheno.perm.mat <- matrix(NA, nrow = num.samples, ncol = n.perm)
	marker.rank.list <- vector(mode = "list", length = dim(pheno)[2])
	marker.rank.mat <- matrix(NA, nrow = n.gene, ncol = n.perm)
	for(i in 1:dim(pheno)[2]){
		marker.rank.list[[i]] <- marker.rank.mat
		}
	names(marker.rank.list) <- colnames(pheno)

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

	one.perm <- function(){
		
	}
	#====================================================


	
	#create a matrix to hold permutation results
  	perm.max <- matrix(NA, ncol = length(pheno[1,]), nrow = n.perm)


      for (j in 1:n.perm) {
      	
      	if(verbose){
 	     	report.progress(current = j, total = n.perm, percent.text = 10, percent.dot = 2)
			}
		
		#shuffle the vector of individuals
        sampled.vector <- sample(1:num.samples)
		gene_perm <- gene[sampled.vector,]
		
		pheno.perm.mat[,j] <- sampled.vector
		
		#loop over traits and do regressions on all markers at once
					
		stat.mat <- matrix(NA, nrow = n.gene, ncol = length(pheno[1,]))
  		for(et in 1:length(pheno[1,])){
  			regress.list <- apply(gene_perm, 2, function(x) lm(pheno[,et]~x))
            stat.mat[,et] <- as.vector(sapply(regress.list, get.stat))
  			}
            
            #find the ranks of the markers across all traits
            max.marker.stat <- apply(stat.mat, 2, rank)
            for(a in 1:dim(max.marker.stat)[2]){
            	marker.rank.list[[a]][,j] <- max.marker.stat[,a]
            	}

                  
      		max.stat <- apply(stat.mat, 2, function(x) max(x, na.rm = TRUE))
            perm.max[j,] <- max.stat
    
        } #end permutations
        
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
