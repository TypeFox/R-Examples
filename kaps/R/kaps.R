##########################################################################
####
#### Multiway-splits adaptive partitioning with distribution free approach 
####
####		Soo-Heang EO and HyungJun CHO 
####
####		Version 1.0.0
####
####		30 Dec, 2013
####
##########################################################################

kaps <- function(formula, data, K = 2:4, mindat, type = c("perm", "NULL"), ...){

	##########################
	##### pre-processing step
	#options(warn = -1)
	if(missing(mindat)) mindat = floor(nrow(data) * 0.05)
	if(!is.Formula(formula)) {
		formula <- as.Formula(formula)
		#class(formula) <- "Formula"
	}

	# for minor options used in kaps algorithm
	# minors = kaps.control()
	minors = kaps.control(...)

	if(any(K == 1)) stop("the number of subgroups (K) have to be greater than subjects.")

	n <- nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n

    if(n == 0L) stop("0 (non-NA) cases.")
	if(length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	if(length(K) > 10) stop("the maximum number of subgroups (K) is too large.")

	call <- match.call()

	######################################################
	###	Finding optimal subgroups
	type <- match.arg(type)

	if(length(K) == 1 & type == "NULL"){
		result = kaps.fit(formula = formula, data = data, K = K, mindat = mindat, minors = minors)
		test.stat2 = adj.test(result)
		result@call <- call
		result@Z <- test.stat2[1] # overall statistic 
		result@X <- test.stat2[2] # test statistic for the worst pair
		return(result)
	} else if(type == "test"){
		cat("Now, finding optimal number of subgroups (K) by test estimates. \n")
		# choose maximal pairwise permutation p-value
		# NEED TO MODIFY
		# lapply(K, kaps.fit, formula = formula, data = data, mindat = mindat, minors = minors)
		test.stat = kaps.test(formula, data, K, minors)
	} else if(type == "perm"){
	
		cat("Now, finding optimal number of subgroups (K) by KAPS-permutation. \n")
		fit = lapply(K, kaps.fit, formula = formula, data = data, mindat = mindat, minors = minors)
		# 1: overall test statistic
		# 2: overall p-value
		# 3: worst pair test statistic
		# 4: worst pair p-value
		test.stat = sapply(fit, kaps.perm, permute = TRUE)
		# 1: overall p-value
		# 2: worst pair p-value
		# choose worst pairwise permutation test

		test.stat2 = sapply(fit, adj.test)
		test.stat2 <- as.matrix(test.stat2)

		index = 1:ncol(test.stat)
		index.pair = test.stat[2,] <= 0.05
		#index.pair <- ifelse(is.na(index.pair), FALSE, index.pair)

		if(all(index.pair == FALSE)){
			# No significant subgroups		
			result <- fit[[1]]
			result@index <- as.integer(0)
			result@groups <- K
			attr(result@groups,"names") <- paste("K<",K,sep="")
		} else{
			index <- index[sum(index.pair)]
			cat("Now, selecting a set of cut-off points...\n")

			result <- fit[[index]]
			result@index <- as.integer(index)
			result@groups <- K
			attr(result@groups,"names") <- paste("K=",K,sep="")
		}

		result@Z <- as.vector(test.stat2[1, ]) #overall statistic 
		result@X <- as.vector(test.stat2[2, ]) #optimal pair p-value
		result@results <- fit # results objects
		result@test.stat <- test.stat # Bonferroni corrected p-values
		result@call <- call
		result@Options <- minors
		return(result)
	
	} else if (type == "NULL"){
		index = 1
	}

	######################################################
	# Obtain log-rank statistics at K
	### parallel computing in order to find optimal k subgroups
	
	cat("Now, selecting a set of cut-off points...\n")
	#fit = kaps.fit(formula = formula, data = data, K = K.sel, mindat = mindat, minors = minors)
	fit = lapply(K, kaps.fit, formula = formula, data = data, mindat= mindat, minors = minors)
	test.stat2 = sapply(fit, adj.test)
	test.stat2 <- as.matrix(test.stat2)

	result <- fit[[index]]
	result@index <- as.integer(index)
	result@groups <- K
	attr(result@groups,"names") <- paste("K=",K,sep="")
	result@Z <- as.vector(test.stat2[1, ]) # overall test statistic 
	result@X <- as.vector(test.stat2[2, ]) # test statistic for the worst pair
	result@results <- fit # result objects
	#if(type != "NULL") {
	#	result@test.stat <- test.stat
	#	result@over.stat.sample <- fit.boot$boot.over.stat
	#	result@pair.stat.sample <- fit.boot$boot.pair.stat
	#}
	result@call <- call
	result@Options <- minors
	return(result)
}

############################################################
# Fit KAPS with test estimates approach 
############################################################
kaps.test <- function(formula, data, K = 2:5, minors = kaps.control()){

	#options(warn = -1)
	
	if(any(K == 1)) stop("the number of subgroups (K) is greater than 1.")

	n = nrow(data) # total number of observations concerning with time 
	rownames(data) <- 1:n

    if(n == 0L) stop("0 (non-NA) cases.")
	if(length(K) < 1) stop("the minimum number of subgroups (K) is greater than 1.")
	#if(length(K) > 10) stop("the maximum number of subgroups (K) is too large.")

	#call = match.call()

	# kaps by test estimates approach
	fold.over.stat = fold.pair.stat = matrix(0, nrow = 1, ncol = length(K))
	result = matrix(0, nrow = length(K), ncol = 4)
	colnames(result) <- c("over_pval", "over_stat", "pair_pval", "pair_stat")
	rownames(result) <- paste("K=",K,sep="")

	## CHECK ME: reduce computing time by parallel computing
	index = sample(1:n, floor(n*0.7))
	learning = data[index,, drop = FALSE]
	test = data[-index,, drop = FALSE]
	mindat = floor(nrow(learning) * 0.05)
	
	fit = lapply(K, kaps.fit, formula = formula, data = learning, mindat = mindat, minors = minors)
	test.stat = sapply(fit, kaps.perm, newdata = test)			
	rownames(test.stat) = c("perm_overall_pval", "perm_min_pair_pval")
	colnames(test.stat) = paste("K=",K,sep="")
	print(round(test.stat,3))

	fold.over.stat[1,] <- test.stat[1,]
	fold.pair.stat[1,] <- test.stat[2,]
	### output
	result[,1] <- fold.over.stat
	#result[,2] <- apply(fold.over.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	result[,3] <- fold.pair.stat
	#result[,4] <- apply(fold.pair.stat, 2, sd, na.rm = TRUE) / sqrt(V)
	#result <- as.data.frame(result)
	return(result)
}
# END @ Nov 12, 2013