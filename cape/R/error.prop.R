error.prop <-
function (data.obj, pairscan.obj, perm = FALSE, verbose = FALSE, n.cores = 2) {

	p = NULL #for appeasing R CMD check

	# require("qpcR")
	# require("corpcor")

	results.obj <- vector(mode = "list", length = 2)
	names(results.obj) <- c("var.to.var.influences.perm")

	if(verbose){
		if(perm){
			cat("\nCalculating error propagation of permuted coefficients.\n")
			}else{
			cat("\nCalculating error propagation of coefficients.\n")	
			}
		}

	#====================================================================================
	#begin internal functions
	#====================================================================================
		get.beta.mat <- function(marker.pair.number){
			#the beta matrix is composed of the coefficients from each pairwise
			#marker model (except for the interaction coefficient)
			#we use all markers in each row and set the non-covariate entries to 0
			num.pheno <- length(scan.two.results)
			beta.mat <- sapply(scan.two.results, function(x) as.numeric(x[[1]][marker.pair.number,(dim(scan.two.results[[1]][[1]])[2]-2):(dim(scan.two.results[[1]][[1]])[2])]))
			rownames(beta.mat) <- c("marker1", "marker2", "interaction")
			return(beta.mat)	
			}


		get.se.mat <- function(marker.pair.number){
			#the beta matrix is composed of the coefficients from each pairwise
			#marker model (except for the interaction coefficient)
			#we use all markers in each row and set the non-covariate entries to 0
			num.pheno <- length(scan.two.results)
			se.mat <- sapply(scan.two.results, function(x) as.numeric(x[[2]][marker.pair.number,(dim(scan.two.results[[1]][[1]])[2]-2):(dim(scan.two.results[[1]][[1]])[2])]))
			rownames(se.mat) <- c("marker1", "marker2", "interaction")
			return(se.mat)	
			}


		
		get.cov.mat <- function(marker.pair.number){
			#the variance-covariance matrix is block diagonal, and 
			#contains three times the number of rows and columns
			#as scanned traits. The blocks contain the variance-
			#covariance matrix from each trait for the marker pair
			cov.mat <- matrix(0, num.pheno*3, num.pheno*3)
			pheno.num <- 1:num.pheno
			start.row.col <- (pheno.num*3)-2
			end.row.col <- (pheno.num*3)
			for(ph in 1:num.pheno){
				cov.mat[start.row.col[ph]:end.row.col[ph], start.row.col[ph]:end.row.col[ph]] <- matrix(scan.two.results[[ph]][[3]][marker.pair.number,], ncol = 3)
				}	
			return(cov.mat)
			}
	
		#check to see if phenotypes or eigentraits were scanned
		pheno.names <- names(pairscan.obj$pairscan.results)	
		pheno.check <- match(pheno.names, colnames(data.obj$pheno))
		if(length(which(!is.na(pheno.check))) == 0){ #if we scanned eigentraits
			num.pheno <- dim(data.obj$ET)[2] #the number of phenotypes
			names.pheno <- colnames(data.obj$ET)
			}else{
			num.pheno <- dim(data.obj$pheno)[2] #the number of phenotypes
			names.pheno <- colnames(data.obj$pheno)
			}				
			
		get.pair.coeffs <- function(marker.pair.number){
			beta.main <- t(get.beta.mat(marker.pair.number)) ### Extract Main effect and interactions
			non.zero <- which(beta.main[1:2,] != 0)
			na.locale <- which(is.na(beta.main))
			if(length(non.zero) > 0 && length(na.locale) == 0){
				beta.se <- t(get.se.mat(marker.pair.number)) ### Extract Main effect and interactions
				beta.cov <- get.cov.mat(marker.pair.number) ### Extract Covars
				inf.coeffs <- calc.delta.errors(markers = marker.mat[marker.pair.number,], beta.m = beta.main, se = beta.se, beta.cov)
				}else{
				inf.coeffs <- NULL	
				}
			return(inf.coeffs)	
			}	
	#====================================================================================
	#end internal functions
	#====================================================================================
		
	### For all marker pairs calculate activity and IC
	
	if(perm){
		if(is.null(pairscan.obj$pairscan.perm)){
			stop("pairscan() with permutations must be run before error.prop()")
			}
		n.perm <- dim(pairscan.obj$pairscan.perm[[1]][[1]])[1]/dim(data.obj$pairscan.results[[1]][[1]])[1]
		marker.mat <-  pairscan.obj$pairscan.perm[[1]][[1]][,1:2]#get all the pairs that were tested in the pair scan
    	scan.two.results <- pairscan.obj$pairscan.perm  #the coefficient matrices from the 2D scan			
		}else{
		if(is.null(pairscan.obj$pairscan.results)){
			stop("pairscan() must be run before error.prop()")
			}
		marker.mat <- pairscan.obj$pairscan.results[[1]][[1]][,1:2] #a matrix listing the names of used marker combinations
	    scan.two.results <- pairscan.obj$pairscan.results  #the coefficient matrices from the 2D scan
		}

    
    colnames(marker.mat) <- c("marker1", "marker2")
	n.pairs <- length(marker.mat[,1]) #number of pairs of genes

	cl <- makeCluster(n.cores)
	registerDoParallel(cl)
	influence.coeffs <- foreach(p = 1:n.pairs, .combine = "rbind") %dopar% {
		get.pair.coeffs(p)
		}				
	stopCluster(cl)
	
	colnames(influence.coeffs) <- c("marker1","marker2","m12","m12.std.dev","m21","m21.std.dev")

	if(perm){
		data.obj$var.to.var.influences.perm <- influence.coeffs
		}else{
		data.obj$var.to.var.influences <- influence.coeffs			
		}		

	if(verbose){
		cat("\n") #make sure the prompt is on the next line at the end of everything
		}
	
	return(data.obj)
}
