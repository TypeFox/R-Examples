pairwise.logrank.test <- function(x, data, formula, rho, adj, splits, shortcut){
	###################################################################
	## Step 1. finding a set of cut-off points 
	## 1-1. Calculate minimum test statistics among possible pairs in the candidate s
	###################################################################
	# treat pre-determined range for predictors
	data$subgroups <- x
	pair <- combnat(unique(x),2)
	#pair <- combn(unique(x),2)
	
	if(shortcut) pair <- pair[, !duplicated(pair[1,])]
	if(is.null(ncol(pair))) pair <- as.matrix(pair)
	
	res <- matrix(NA, nrow = ncol(pair), ncol=3)
	for(i in 1:ncol(pair)){
		data.tmp <- data[data$subgroups %in% pair[,i],]
		if(splits == "logrank"){
			tmp <- try(survdiff(formula = formula, data = data.tmp, rho = rho), silent = TRUE)
			if(class(tmp) == "try-error"){
				res[i,1] <- 0 # pair-wise test statistic
				res[i,2] <- 1 # pair-wise p-value
				res[i,3] <- as.numeric(paste(pair[,i], collapse ="")) # pairs for candidates
			} else {
				res[i,1] <- tmp$chisq # pair-wise test statistic
				res[i,2] <- 1 - pchisq(tmp$chisq, 1) # pair-wise p-value
				res[i,3] <- as.numeric(paste(pair[,i], collapse ="")) # pairs for candidates				
			}
		} #else if(splits == "exact"){
			#class(formula) <- "formula"
			#tmp = surv_test(formula = formula, data = data.tmp, distribution = exact())
			#tmp <- maxstat.test(formula= formula, data = data.tmp, smethod="LogRank",
			#	pmethod="exactGauss",minprop = 0.01, maxprop = .99)
			#res[i,1] <- tmp$statistic # pair-wise test statistic
			#res[i,2] <- tmp$p.value # pair-wise p-value
			#res[i,1] <- tmp@statistic@teststatistic
			#res[i,2] <- 1 - pchisq(res[i,1], tmp@statistic@df)
			#if(is.na(res[i,2])) res[i,2] <- 0
			#res[i,3] <- as.numeric(paste(pair[,i], collapse ="")) # pairs for candidates
		#}
	}
	# Adjsut P-values for Multiple Comparisons 
	if(adj != "none") res[,2] <- p.adjust(res[,2], method = adj)
	#if(cre == "pvalue") res <- res[which(res[,2] == max(res[,2])),,drop = FALSE]
	#else if(cre == "chisq") res <- res[which(res[,1] == min(res[,1])),, drop = FALSE]
	index <- which(res[,1] == min(res[,1], na.rm = TRUE))
	res <- res[index,,drop = FALSE]
	if(nrow(res) >=2) res = res[1,, drop = FALSE] 
	return(res)
}
 
group.sel <- function(x.vec, pt, K, mindat, data, f, minors){
	###################################################################
	## Step 1. finding a set of cut-off points 
	## 1-1. Calculate worst pair test statisic among possible pairs in the candidate s
	## 1-2. Select a representative candidate pair with the largest test statistics
	###################################################################
	### find their groups automatically
	if(length(unique(x.vec)) < K) stop("The # of unique x must be larger than # of groups.")
	candid.pt <- combnat(pt, (K-1))
	#candid.pt <- combn(pt, (K-1))
	nr <- length(x.vec)
	cfun <- function(candid, nr, x.vec, mindat,formula){
		nc <- length(candid)
		gClass <- matrix(NA, ncol = nc, nrow = nr)
		gClass <- sapply(candid, function(x, x.vec) x.vec > x, x.vec = x.vec)
		if(is.vector(gClass)) gClass <- t(gClass)
		where <- apply(gClass, 1, sum)
		where <- where + 1
		if(length(unique(where)) != K) return(c(NA,NA,NA))
		if(all(table(where) > mindat)){
			test.tmp <- pairwise.logrank.test(where, data, formula, 
				rho = minors@rho, 
				adj = minors@p.adjust.methods, 
				splits = minors@splits, 
				shortcut = minors@shortcut)
			return(test.tmp)
		} else {
			return(c(NA,NA,NA))
		}
	}
	result <- apply(candid.pt, 2, cfun, nr = nr, x.vec = x.vec, mindat = mindat, formula = f) 
	if(all(is.na(result))) {
		cat("You have to modify the arguments, K (less than",K,") or mindat (greater than",mindat,"). \n ")
		stop("This parameters does not meet the minimum sample rule in the subgroup.")
	}

	index <- which(result[1,] == max(result[1,], na.rm = TRUE))
	if(length(index) >= 2) index <- index[1]

	##########################
	##allocate subgroup vector
	gClass <- matrix(NA, ncol = K, nrow = nr)
	gClass <- sapply(candid.pt[,index], function(x,y) y > x, y = x.vec)
	where <- apply(gClass, 1, sum) + 1
	gc(reset = TRUE)
	return(list(index = index, test = result, where = where))
}
