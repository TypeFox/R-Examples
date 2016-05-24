calc.p <-
function(data.obj, pairscan.obj, pval.correction = c("holm", "fdr", "lfdr", "none"), n.cores = 2) {
	
	p = NULL #for appeasing R CMD check
	
	choice <- 0
	if(!is.null(data.obj$var.to.var.p.val)){
		message("\nIt appears that p-values have already been calculated for this data object.\n", sep = "")

		while(choice != "a" && choice != "b" && choice != "c"){
			choice <- readline(prompt = "Would you like to:\n\t(a) recalculate all p-values?\n\t(b) change the multiple testing correction?\n\t(c) cancel?\n")
			}
		if(choice == "c"){
			return(data.obj)
			}
		}
	
	plot.null.transform = TRUE
	
	if(length(grep("h", pval.correction) > 0)){
		pval.correction <- "holm"
		}
		
	if(pval.correction != "holm" && pval.correction != "fdr" && pval.correction != "lfdr" && pval.correction != "none"){
		stop("pval.correction must be one of the following: 'holm', 'fdr', 'lfdr', 'none'")
		}
		
	
	if(choice == "a" || choice == 0){
	influences.org <- data.obj$var.to.var.influences
	influences.perm <- data.obj$var.to.var.influences.perm
		
	if(is.null(influences.org)){
		stop("error.prop() with perm = FALSE must be run before running calc.p()")
		}
		
	if(is.null(influences.perm)){
		stop("error.prop() with perm = TRUE must be run before running calc.p()")
		}


	n.gene <- dim(data.obj$geno.for.pairscan)[2] #get the number of genes used in the pair scan
	n.pairs <- dim(influences.org)[1] #the number of pairs with successful calculations
    
    	
    marker.mat <- influences.org[,1:2] #a matrix listing the names of all marker combinations
    colnames(marker.mat) <- c("marker1", "marker2")


	#### Combinine across permutations#####
	#get the t statistics for all permutations
	mat12.perm <- as.numeric(influences.perm[,3]) / as.numeric(influences.perm[,4])
	mat21.perm <- as.numeric(influences.perm[,5]) / as.numeric(influences.perm[,6])
	mat12.mat21.perm <- c(mat12.perm, mat21.perm)

	mat12 <- as.numeric(influences.org[,3]) / as.numeric(influences.org[,4])
	mat21 <- as.numeric(influences.org[,5]) / as.numeric(influences.org[,6])


	adj.null <- mat12.mat21.perm
	adj.m12 <- mat12
	adj.m21 <- mat21


	#changed calculation of p value to account for the asymmetric m12/m21 distribution
	#I now calculate the p value based on above and below the median m12/m21
	#separately
	get.emp.p <- function(num.pair){
		emp.vals <- c(adj.m12[num.pair], adj.m21[num.pair])
		emp.p <- rep(NA, 2)

		for(e in 1:length(emp.vals)){
			if(!is.na(emp.vals)[e]){
				if(emp.vals[e] < median(adj.null)){
					emp.p[e] <- length(which(adj.null <= emp.vals[e]))/length(which(adj.null <= median(adj.null)))
					}else{
					emp.p[e] <- length(which(adj.null >= emp.vals[e]))/length(which(adj.null >= median(adj.null)))
					}
				}else{
				emp.p[e] <- NA	
				}
			}
		return(emp.p)
		}

	cl <- makeCluster(n.cores)
	registerDoParallel(cl)
	all.emp.p <- foreach(p = 1:n.pairs, .combine = "rbind") %dopar% {
		get.emp.p(p)
		}
	stopCluster(cl)
	
	# all.emp.p <- t(apply(matrix(1:n.pairs, ncol = 1), 1, function(x) get.emp.p(x)))

	m12 <- matrix(c(marker.mat[,2],marker.mat[,1],as.numeric(as.matrix(influences.org[,3])),as.numeric(as.matrix(influences.org[,4])),(abs(as.numeric(influences.org[,3])) / as.numeric(influences.org[,4])),all.emp.p[,1]), ncol = 6)	
	colnames(m12) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")	

	m21 <- matrix(c(marker.mat[,1],marker.mat[,2],as.numeric(as.matrix(influences.org[,5])),as.numeric(as.matrix(influences.org[,6])),(abs(as.numeric(influences.org[,5])) / as.numeric(influences.org[,6])),all.emp.p[,2]), ncol = 6)
	colnames(m21) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")
	} #end if choice == "a"
	
	
	#adjust the p values
	if(choice == "b"){
		final.table <- data.obj$var.to.var.p.val
		final.table <- final.table[,-7]
		}else{
		final.table <- rbind(m12, m21)
		}
		
	if(pval.correction == "none"){
		p.adjusted <- as.numeric(final.table[,"P_empirical"])
		final.table <- cbind(final.table, p.adjusted)
		}
	if(pval.correction == "holm"){
		p.adjusted <- p.adjust(as.numeric(final.table[,"P_empirical"]), method = "holm")
		final.table <- cbind(final.table, p.adjusted)
		}
		
	if(pval.correction == "fdr" || pval.correction == "lfdr"){
		pvals <- as.numeric(final.table[,"P_empirical"])
		not.na.locale <- which(!is.na(pvals))
		fdr.out <- fdrtool(as.numeric(final.table[not.na.locale,"P_empirical"]), statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "fndr")
		if(pval.correction == "lfdr"){
			lfdr <- rep(NA, dim(final.table)[1])
			lfdr[not.na.locale] <- fdr.out$lfdr
			final.table <- cbind(final.table, lfdr)
			}else{
			qval <- rep(NA, dim(final.table)[1])
			qval[not.na.locale] <- fdr.out$qval
			final.table <- cbind(final.table, qval)	
			}
		}


	final.table <- final.table[order(final.table[,"|Effect|/SE"], decreasing = TRUE),]

	data.obj$var.to.var.p.val <- final.table
	
	return(data.obj)
}
