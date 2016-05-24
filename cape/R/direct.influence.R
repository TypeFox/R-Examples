direct.influence <-
function(data.obj, pairscan.obj, transform.to.phenospace = TRUE, pval.correction = c("holm", "fdr", "lfdr"), verbose = FALSE, save.permutations = FALSE) {


	#This object keeps large intermediate data
	intermediate.data <- vector(mode = "list")

		
	if(length(grep("h", pval.correction) > 0)){
		pval.correction = "holm"
		}

	if(pval.correction != "holm" && pval.correction != "fdr" && pval.correction != "lfdr" && pval.correction != "none"){
		stop("Method must be one of the following: 'holm', 'fdr', 'lfdr', 'none'")
		}
		

	data.obj$transform.to.phenospace <- transform.to.phenospace


	#calculate the direct influences for either the actual
	#tests or the permutations. Return a list with one element
	#for each phenotype with columns:
	#marker1, marker2, marker1.influence.coef, marker2.influence.coef, marker1.se, marker2.se
	
	dir.inf <- function(data.obj, perm){
		geno <- data.obj$geno.for.pairscan 
		n.gene <- dim(geno)[2]
		if(perm){
			scan.two.results <- pairscan.obj$pairscan.perm		
			intermediate.data$pairscan.results <- scan.two.results
			pairscan.obj$pairscan.perm <- NULL
			}else{
			scan.two.results <- pairscan.obj$pairscan.results
			}
			
		n.pairs <- length(scan.two.results[[1]][[1]][,1]) 
		marker.mat <- scan.two.results[[1]][[1]][,1:2] 
		colnames(marker.mat) <- c("marker1", "marker2")
		orig.pheno <- data.obj$pheno
	
		
		#if transform.to.phenospace is FALSE, figure out
		#if we are calculating the direct influences on
		#phenotypes or eigentraits
	
		if(!transform.to.phenospace){
			pheno.names <- names(pairscan.obj$pairscan.results)	
			pheno.check <- match(pheno.names, colnames(data.obj$pheno))
			if(length(which(!is.na(pheno.check))) == 0){ #if we scanned eigentraits
				ET <- data.obj$ET					
				}else{
				ET <- data.obj$pheno	
				}					
			}else{
			ET <- data.obj$ET
			right.sing.vals <- data.obj$right.singular.vectors 
			diag.mat <- round(t(ET)%*%orig.pheno%*%right.sing.vals, 2) #get the singular value matrix
			pheno.names <- colnames(data.obj$pheno)			
			}			

			
		num.ET <- dim(ET)[2] 
		num.pheno <- length(pheno.names)
		coeff.names <- c("marker1", "marker2")
		
		#=========================================================
		# preallocate matrices to hold the statistics for each pair
		# The columns are for marker1, marker2, marker1.beta, 
		# marker2.beta, marker1.se and marker2.se (6 columns)
		# there is one of these for each phenotype
		#=========================================================
		stats.mat <- matrix(NA, ncol = 6, nrow = n.pairs)
		colnames(stats.mat) <- c("marker1", "marker2", "marker1.coef", "marker2.coef", "marker1.se", "marker2.se")
		stats.list <- vector(mode = "list", length = num.pheno)


		names(stats.list) <- pheno.names
		for(n in 1:num.pheno){
			stats.list[[n]] <- stats.mat
			}


		#This function grabs either the beta matrix or the se matrix
		#The result element determines which type of matrix will be 
		#returned: beta, result.element = 1, se, result.element = 2
		
		get.beta.prime <- function(scan.two.results, marker.pair.number){

			#the beta matrix is composed of the coefficients from each pairwise
			#marker model (except for the interaction coefficient)
			#we use all markers in each row and set the non-covariate entries to 0
			num.pheno <- length(scan.two.results)
			# beta.mat <- matrix(NA, nrow = (dim(scan.two.results[[1]][[1]])[2]-3), ncol = num.pheno)
			beta.mat <- matrix(NA, nrow = 2, ncol = num.pheno)
			for(ph in 1:num.pheno){
				beta.mat[,ph] <- as.numeric(scan.two.results[[ph]][[1]][marker.pair.number,3:(dim(scan.two.results[[ph]][[1]])[2]-1)])
				}				
			colnames(beta.mat) <-colnames(ET)
			rownames(beta.mat) <- colnames(scan.two.results[[ph]][[1]])[3:(dim(scan.two.results[[ph]][[1]])[2]-1)]
			

			#calculate the full beta prime matrix
			#if we are converting to phenotype space
			if(transform.to.phenospace){
				beta.prime <- beta.mat%*%diag.mat%*%t(right.sing.vals)
				#get only the stats for marker1 and marker2
				trunc.beta.prime <- matrix(beta.prime[(length(beta.prime[,1])-1):length(beta.prime[,1]),], nrow = 2)
				colnames(trunc.beta.prime) <- colnames(orig.pheno)
				}else{
				beta.prime <- beta.mat #otherwise, just use the straight beta matrix
				trunc.beta.prime <- matrix(beta.prime[(length(beta.prime[,1])-1):length(beta.prime[,1]),], nrow = 2)
				colnames(trunc.beta.prime) <- colnames(ET)
				}

			just.marker.beta <- vector(mode = "list", length = dim(trunc.beta.prime)[2])
			for(i in 1:length(just.marker.beta)){
				just.marker.beta[[i]] <- trunc.beta.prime[,i]
				}
			names(just.marker.beta) <- colnames(trunc.beta.prime)
			return(just.marker.beta)
	
			}


		get.se.prime <- function(scan.two.results, marker.pair.number){
			#also get the se prime matrix
			num.pheno <- length(scan.two.results)
			se.mat <- matrix(NA, nrow = (dim(scan.two.results[[1]][[1]])[2]-3), ncol = num.pheno)
			for(ph in 1:num.pheno){
				se.mat[,ph] <- as.numeric(scan.two.results[[ph]][[2]][marker.pair.number,3:(dim(scan.two.results[[ph]][[2]])[2]-1)])
				}				
			colnames(se.mat) <-colnames(ET)
			rownames(se.mat) <- colnames(scan.two.results[[ph]][[2]])[3:(dim(scan.two.results[[ph]][[2]])[2]-1)]
			
	
			#Calculate a beta se prime matrix from each
			#beta matrix by multiplying by the beta se, squared diagonal singular value 
			#matrix, and the square of the transpose of the right singular values.	

			sqrd.se.mat <- se.mat ^ 2
			if(transform.to.phenospace){
				temp.mat <- diag.mat %*% t(right.sing.vals)
				sqrd.temp.mat <- (temp.mat ^  2)
				se.prime <- sqrd.se.mat %*% sqrd.temp.mat
				colnames(se.prime) <- colnames(orig.pheno)
				}else{
					se.prime <- sqrd.se.mat
					colnames(se.prime) <- colnames(ET)
					}

			se.prime <- sqrt(se.prime)

			just.marker.se <- vector(mode = "list", length = dim(se.prime)[2])
			for(i in 1:length(just.marker.se)){
				just.marker.se[[i]] <- se.prime[((dim(se.prime)[1]-1):dim(se.prime)[1]),i]
				}
			names(just.marker.se) <- colnames(se.prime)
			return(just.marker.se)
	
			}
		
	

		#get the beta matrix and se matrix for each pair of markers
		all.beta.prime <- apply(matrix(1:dim(marker.mat)[1], ncol = 1), 1, function(x) get.beta.prime(scan.two.results, x))
					
		#add these into the final stats.matrix
		for(i in 1:length(stats.list)){
			stats.list[[i]][,1:2] <- marker.mat
			stats.list[[i]][,3:4] <- t(sapply(all.beta.prime, function(x) x[[i]]))
			}
		all.beta.prime <- NULL
		
		#also get the se matrices
		all.se.prime <- apply(matrix(1:dim(marker.mat)[1], ncol = 1), 1, function(x) get.se.prime(scan.two.results, x))
					
		#add these into the final stats.matrix
		for(i in 1:length(stats.list)){
			stats.list[[i]][,5:6] <- t(sapply(all.se.prime, function(x) x[[i]]))
			}
		all.se.prime = NULL
	
			
			return(stats.list)
			}


		#================================================================
		# calculate the direct influence of the variants and permutations
		#================================================================
			
		if(verbose){cat("calculating direct influence of variants...\n")}
		exp.stats <- dir.inf(data.obj, perm = FALSE)
		if(verbose){cat("calculating direct influence of permutations...\n")}
		perm.stats <- dir.inf(data.obj, perm = TRUE)
		

		intermediate.data$var.to.pheno.influence <- exp.stats
		intermediate.data$var.to.pheno.influence.perm <- perm.stats

		#================================================================
		#================================================================	



	direct.influence.stat <- function(pair.stats, perm){
			

		#separate out the markers from each pair
		#for each marker in each context, give it an 
		#influence coefficient, influence se, t stat
		#and |t stat|

		markers <- sort(unique(as.numeric(c(pair.stats[[1]][,1], pair.stats[[1]][,2]))))
		pheno.names <- names(pair.stats)
		
		
		marker.incidence <- apply(matrix(markers, ncol = 1), 1, function(x) length(which(pair.stats[[1]][,1:2] == x)))
		#preallocate a matrix that will hold the statistics for each marker
		#in each pair context. The columns are marker, coef, se, t.stat
		#|t.stat|, emp.p (6 columns)
		if(perm){
			stats.mat <- matrix(NA, nrow = sum(marker.incidence), ncol = 5)
			colnames(stats.mat) <- c("marker", "coef", "se", "t.stat", "|t.stat|")
			}else{
			stats.mat <- matrix(NA, nrow = sum(marker.incidence), ncol = 6)
			colnames(stats.mat) <- c("marker", "coef", "se", "t.stat", "|t.stat|", "emp.p")			
			}
		
		ind.stats.list <- vector(mode = "list", length = length(pheno.names))
		names(ind.stats.list) <- pheno.names
		for(i in 1:length(ind.stats.list)){
			ind.stats.list[[i]] <- stats.mat
			}
		
		#for each phenotype, go through the marker names and
		#collect the statistics from var.to.pheno.influence 
		#for when the marker was in position 1 and in position 2
		for(ph in 1:length(pheno.names)){
			
			#go through each marker and take out its statistics
			for(m in markers){

				#find out where the marker was in position 1
				m1.locale <- which(pair.stats[[ph]][,1] == m)
				
				if(length(m1.locale) > 0){
					beta.coef1 <- as.numeric(pair.stats[[ph]][m1.locale, "marker1.coef"])
					se1 <- as.numeric(pair.stats[[ph]][m1.locale, "marker1.se"])
					t.stat1 <- beta.coef1/se1
					stat.section1 <- matrix(c(rep(as.numeric(m), length(beta.coef1)), beta.coef1, se1, t.stat1, abs(t.stat1)), ncol = 5)

					start.row <- min(which(is.na(ind.stats.list[[ph]][,1])))
					ind.stats.list[[ph]][start.row:(start.row + dim(stat.section1)[1] - 1),1:dim(stat.section1)[2]] <- stat.section1
					}

				#find out where the marker was in position 2
				m2.locale <- which(pair.stats[[ph]][,2] == m)

				if(length(m2.locale) > 0){
					beta.coef2 <- as.numeric(pair.stats[[ph]][m2.locale, "marker2.coef"])
					se2 <- as.numeric(pair.stats[[ph]][m2.locale, "marker2.se"])
					t.stat2 <- beta.coef2/se2
					stat.section2 <- matrix(c(rep(as.numeric(m), length(beta.coef2)), beta.coef2, se2, t.stat2, abs(t.stat2)), ncol = 5)
				
					start.row <- min(which(is.na(ind.stats.list[[ph]][,1])))
					ind.stats.list[[ph]][start.row:(start.row + dim(stat.section2)[1] - 1),1:dim(stat.section2)[2]] <- stat.section2
					}
	
				}

			}
			
			return(ind.stats.list)
	
			}


		#================================================================
		# tabulate influences of individual markers		
		#================================================================
		if(verbose){cat("calculating individual marker influences...\n")}
		exp.stats.per.marker <- direct.influence.stat(pair.stats = exp.stats, perm = FALSE)
		if(verbose){cat("calculating individual marker influences in permutations...\n")}
		perm.stats.per.marker <- direct.influence.stat(perm.stats, perm = TRUE)
		
		#================================================================
	
	direct.influence.epcal<- function(exp.stats.per.marker, perm.stats.per.marker){
	
		stat <- exp.stats.per.marker
		perm.test.stat <- perm.stats.per.marker
		
		get.p <- function(val, dist){
			pval <- length(which(dist >= val))/length(dist)
			return(pval)
			}

		#caclualte an empirical p value for each entry in stat
		for(ph in 1:length(stat)){
			comp.dist <- as.numeric(perm.test.stat[[ph]][,"|t.stat|"]) #compare each calculated value to the permuted distribution
			stat[[ph]][,"emp.p"] <- apply(matrix(as.numeric(stat[[ph]][,"|t.stat|"]), ncol = 1), 1, function(x) get.p(x, comp.dist))
			}
		
		return(stat)
		}


		#================================================================
		# calculate empirical p values
		#================================================================
		if(verbose){cat("calculating p values...\n")}
		emp.p.table <- direct.influence.epcal(exp.stats.per.marker, perm.stats.per.marker)
		#replace the direct influence tables with the tables with the added empirical p values
		intermediate.data$var.to.pheno.test.stat <- emp.p.table
		intermediate.data$var.to.pheno.test.stat.perm <- perm.stats.per.marker
		
		#================================================================
		#================================================================
		

	direct.influence.max <- function(emp.p.table){
		
		stat <- emp.p.table

		markers <- unique(stat[[1]][,1])
		max.stat.list <- vector(length(stat), mode = "list")
		names(max.stat.list) <- names(stat)
		
		#go through each of the phenotypes and find the maximum 
		#influence of each marker
		for(ph in 1:length(stat)){
			max.stat.table <- NULL
			for(m in markers){
				marker.locale <- which(stat[[ph]][,1] == m)
				temp.table <- matrix(stat[[ph]][marker.locale,], nrow = length(marker.locale))
				colnames(temp.table) <- colnames(stat[[ph]])
				max.stat.locale <- which(as.numeric(temp.table[,"|t.stat|"]) == max(as.numeric(temp.table[,"|t.stat|"]), na.rm = TRUE))
				if(pval.correction == "holm"){
					adj.p <- p.adjust(temp.table[,6], method = "holm")
					adj.p.name <- "p.adjusted"
					}
				if(pval.correction == "fdr"){
					fdr.out <- suppressWarnings(fdrtool(temp.table[,6], statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "pct0"))
					adj.p <- fdr.out$qval
					adj.p.name <- "qval"
					}
				if(pval.correction == "lfdr"){
					fdr.out <- suppressWarnings(fdrtool(temp.table[,6], statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "pct0"))
					adj.p <- fdr.out$lfdr
					adj.p.name <- "lfdr"
					}
				max.stats <- matrix(c(temp.table[max.stat.locale,], rep(min(adj.p), length(max.stat.locale))), nrow = length(max.stat.locale))
				max.stat.table <- rbind(max.stat.table, max.stats[1,]) #if there is more than one max, just take one
				colnames(max.stat.table) <- c(colnames(stat[[1]]), adj.p.name)
				}	
			ordered.table <- max.stat.table[order(max.stat.table[,5], decreasing = TRUE),]
			max.stat.list[[ph]] <- ordered.table
			}
	
		return(max.stat.list)
		}


	if(save.permutations){
		saveRDS(intermediate.data, "permutation.data.RData")
		}	

	if(verbose){cat("calculating maximum influence of each marker...\n")}
	max.stat.list <- direct.influence.max(emp.p.table)
	data.obj$max.var.to.pheno.influence <- max.stat.list
	
	
	return(data.obj)
}
