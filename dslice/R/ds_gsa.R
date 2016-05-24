ds_gsa <- function(expdat, geneset, label, generank, ..., lambda = 1, bycol = FALSE,
                   minsize = 15, maxsize = 500, randseed = 11235, rounds = 1000)
{
	if(is.character(expdat)){
		expdat <- load_gct(expdat)
	}
	if(is.character(geneset)){
		geneset <- load_gmt(geneset)
	}
	if(is.character(label)){
		label <- load_cls(label)
	}
	n_gene <- nrow(expdat)
	n_sample <- ncol(expdat)
	n_sets <- length(geneset$set_name)
	level <- max(label$value) + 1
	if(n_sample != length(label$value)){
		stop("label do not match expression data")
	}
	if(!is.numeric(generank)){
		if(is.character(generank)){
			rankmethod <- get(generank, mode = "function", envir = parent.frame())
		}
		if(!is.function(rankmethod)){
			stop("'generank' must be a integer vector contains valid rank or a function or a string naming a valid function.")
		}
		generank <- rankmethod(expdat, label$value, ...)
	}
	if(length(generank) != n_gene){
		stop("rank list does not match expression data")
	}
	genename <- rownames(expdat)
	sortname <- genename[generank]
	affmat <- matrix(rep(0, n_sets * n_gene), nrow = n_sets, ncol = n_gene)
	used_set_num <- 1
	used_set_name <- rep("null", n_sets)
	used_set_syb_num <- rep(0, n_sets)
	used_set_symbol <- NULL
	symbol_list <- geneset$gene_symbol
	for(i in 1:n_sets){
		set_gene <- symbol_list[[i]]
		exist_set_gene <- intersect(set_gene, genename)
		exist_set_size <- length(exist_set_gene)
		if((exist_set_size < minsize) || (exist_set_size > maxsize)){
			next
		}
		used_set_name[used_set_num] <- geneset$set_name[i]
		used_set_syb_num[used_set_num] <- exist_set_size
		used_set_symbol[[used_set_num]] <- exist_set_gene
		affmat[used_set_num, match(exist_set_gene, genename)] <- 1L
		used_set_num <- used_set_num + 1
	}
	used_set_num <- used_set_num - 1
	affmat <- affmat[1:used_set_num, ]
	used_set_name <- used_set_name[1:used_set_num]
	used_set_syb_num <- used_set_syb_num[1:used_set_num]
	# slicing
	pflag <- rep(0, used_set_num)  #  indicating which one is picked out (DS-value > 0)
	ds_val <- rep(0, used_set_num)  #  DS-value, nonnegative
	sortidx <- match(1:n_gene, generank)
	sortmat <- affmat[, sortidx]
	slices <- NULL
	for(i in 1:used_set_num){
	current_set <- as.integer(sortmat[i, ])
	current_res <- ds_k(current_set, level, lambda, slice=TRUE)
	if(current_res$dsval > 1e-6){
		pflag[i] <- 1
		ds_val[i] <- current_res$dsval
	}
	# record slices
	tmpslice <- current_res$slices
	colnames(tmpslice) <- c(label$pheotype, "total")
	slices[[i]] <- tmpslice
	}
	nullval <- matrix(1, nrow = used_set_num, ncol = rounds)
	set.seed(seed = randseed, kind = NULL)
	for(j in 1:rounds){
		if(bycol){
			permrank <- rankmethod(expdat, sample(label$value), ...)
		}else{
			permrank <- sample(generank)
		}
		permidx <- match(1:n_gene, permrank)
		permmat <- affmat[, permidx]
		for(i in 1:used_set_num){
			current_set <- as.integer(permmat[i, ])
			nullval[i, j] <- ds_k(current_set, level, lambda)
		}
	}
	pvalues <- vector(length = used_set_num, mode = "numeric")
	qvalues <- rep(0, used_set_num)

	for(i in 1:used_set_num){
		current_val <- ds_val[i]
		if(pflag[i] == 1){
			pvalues[i] <- length(which(nullval[i, ] > current_val))
		}else{
			pvalues[i] <- rounds
		}
	}
	pvalues <- pvalues / rounds
	qvalues <- pvalues * used_set_num / rank(pvalues)
	qvalues <- pmin(qvalues, 1)
	GSA <- list(
		set_name = used_set_name, 
		set_size = used_set_syb_num,
		DS_value = ds_val, 
		pvalue = pvalues, 
		FDR = qvalues, 
		slices = slices
	)
	return(GSA)
}
