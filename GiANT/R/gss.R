##########################################################################
#gene set statistics
#
# x:				One value per gene gained from gls/transformation.
# geneSetIndices:	The indices of the genes in the gene set in x.
##########################################################################
gss.mean <- function(x, geneSetIndices){
	return(mean(x[geneSetIndices]))
}

gss.sum <- function(x, geneSetIndices){
	return(sum(x[geneSetIndices]))
}

gss.wilcoxonRankTest <- function(x, geneSetIndices){
	tmp <- wilcox.test(x[geneSetIndices])
	return(tmp$p.value)
}

gss.maxmean <- function(x, geneSetIndices){
	tmp.x <- x[geneSetIndices]

	pos <- (sum(tmp.x>=0)*mean(tmp.x[tmp.x>=0])) / length(geneSetIndices)
	neg <- (-sum(tmp.x<0)*mean(tmp.x[tmp.x<0])) / length(geneSetIndices)

	return(max(pos,neg))
}

gss.median <- function(x, geneSetIndices){
	return(median(x[geneSetIndices]))
}

gss.enrichmentScore <- function(x, geneSetIndices, p = 1){
	
	if(sum((x == Inf) | (x == -Inf)) != 0){
		warning("Gene set statistic: enrichmentScore: Inf values found.\nInf is set to max(x), -Inf is set to min(x).")

		inf_ind <- which(x == Inf)
		minf_ind <- which(x == -Inf)

		x[x == Inf] <- max(x[-c(inf_ind,minf_ind)])
		x[x == -Inf] <- min(x[-c(inf_ind,minf_ind)])
	}

	indices <- order(x, decreasing=TRUE)

	NR <- sum(abs(x[geneSetIndices])^p)

	hits <- rep(0,length(x))
	hits[geneSetIndices] <- (x[geneSetIndices]^p) / NR 
	hits <- hits[indices]

	misses <- rep(1/(length(x)-length(geneSetIndices)),length(x))
	misses[geneSetIndices] <- 0
	misses <- misses[indices]

	ES <- cumsum(hits - misses)
	maxES <- ES[which.max(abs(ES))]
	
	return(maxES)
}

gss.fisherExactTest <- function(x, geneSetIndices){

	nAllGenes <- length(x)
	coreSet <- which(x != 0)

	##########################################
	# 2 by 2 table
	##########################################
	tp <- sum(is.element(geneSetIndices, coreSet))
	# if tp == 0 the p-Value will be 1. No further calculations are necessary
	if(tp == 0){
		return(list(
			pValue = 1,
			coreSet = coreSet,
			intersectGeneSetCoreSet = NULL,
			nAllGenes = nAllGenes,
			res.all = NULL))
	}
	##########################################
	tn <- sum(!is.element(1:nAllGenes, c(coreSet, geneSetIndices)))
	fp <- length(geneSetIndices) - tp
	fn <- sum(!is.element(coreSet, geneSetIndices))

	table_2by2 <- matrix(c(tp, fp, fn, tn),2,2)

	##########################################
	# R fisher test
	##########################################
	f.res <- fisher.test(table_2by2, alternative = "greater")

	return(list(
		pValue = f.res$p.value,
		coreSet = names(x)[coreSet],
		intersectGeneSetCoreSet = names(x)[intersect(geneSetIndices, coreSet)],
		nAllGenes = nAllGenes,
		res.all = f.res))
}


#### The Gene Set Z-score
## Toeroenen et al 2009:
# Robust extraction of functional signals from gene set analysis
# using a generalized threshold free scoring function
gss.gsz <- function(x,
	geneSetIndices,
	w1 = 0.2,
	w2 = 0.5,
	preVar = 0,
	varConstant = 10){

	tmpGeneSets <- rep(0,length(x))
	tmpGeneSets[geneSetIndices] <- 1
	#duplicate gene set because of missing "drop=FALSE" in mGSZ
	tmpGeneSets <- matrix(tmpGeneSets, ncol = 1, nrow = length(tmpGeneSets))

	pos.scores <- GSZ.test.score(
		expr.data=x,
		GO.pointers = tmpGeneSets,
		wgt1 = w1,
		wgt2 = w2,
		pre.var = preVar,
		prior.cl.sz = varConstant)

#	print(pos.scores)

#	pos.mGSZ.scores <- pos.scores$gene.set.scores$mGSZ.scores

	return(pos.scores$result)
}


#############################################################################################
# Adapted from http://ekhidna.biocenter.helsinki.fi/users/petri/public/GSZ/GSZ_code_test.R #
#############################################################################################
GSZ.test.score <- function(expr.data, GO.pointers, wgt1 = 0.2, wgt2 = 0.5, prior.cl.sz = 10, pre.var = 0){

	num_genes <- length(expr.data)

	ord_out <- order(expr.data, decreasing= TRUE)
	expr.data <- expr.data[ord_out]
	go.sz <- dim(GO.pointers)
	cols <- go.sz[2]
	GO.pointers <- GO.pointers[ord_out,,drop=FALSE]

	expr.data.ud <- expr.data[num_genes:1]	# expression values turned up-side down
											# This does the analysis of the lower end

	# Define mean values and squared gene expression mean values
	divider <- c(1:num_genes)
	mean_table    <- cumsum(expr.data)/divider
	mean_table_ud <- cumsum(expr.data.ud)/divider
	mean_table_sq <- mean_table^2
	mean_table_ud_sq <- mean_table_ud^2

	# Define the variance values    
	#pre_var <- 0	# This is the uncertainty of the observations. A value often omitted from variance
					# Omitted also here (value = 0), but included for potential usage

	var_table    <- cumsum(expr.data^2)/divider - (mean_table)^2 + pre.var
	var_table_ud <- cumsum(expr.data.ud^2)/divider - (mean_table_ud)^2 + pre.var 

	max_val <- c(1:num_genes)

	# Define prior variances... 
	tmp1 <- GSZ_calc_prior(var_table, mean_table_sq, num_genes, prior.cl.sz)
	tmp2 <- GSZ_calc_prior(var_table_ud, mean_table_ud_sq, num_genes, prior.cl.sz)
	var_const1 <- wgt1*tmp1
	var_const2 <- wgt1*tmp2

	A <- matrix(0, num_genes, cols)
	B <- A

	for(k in 1:cols){

		tmp <- c(1:num_genes)
		po1 <- which(GO.pointers[,k] == 1)
		po0 <- tmp[-po1]

		if(length(po1) > 0 && length(po0) > 0) { 
			hyge.stat <- hyge.mean.var(num_genes, length(po1))
			Z_mean1 <- mean_table*(2*hyge.stat$mean - max_val)
			Z_mean2 <- mean_table_ud*(2*hyge.stat$mean - max_val)
			prob_sum <- count.prob.sum(num_genes, hyge.stat$mean, hyge.stat$var)
			Z_var1 <- 4*(var_table*prob_sum    + mean_table_sq*hyge.stat$var)
			Z_var2 <- 4*(var_table_ud*prob_sum + mean_table_ud_sq*hyge.stat$var)
			tmp1 <- expr.data
			tmp1[po0] <- 0
			tmp0 <- expr.data
			tmp0[po1] <- 0
			tulos1 <- cumsum(tmp1) - cumsum(tmp0) - Z_mean1
			tulos2 <- cumsum(tmp1[num_genes:1]) - cumsum(tmp0[num_genes:1]) - Z_mean2

			A[,k] <- tulos1/(Z_var1 + wgt2*median(Z_var1) + var_const1)^0.5
			B[,k] <- tulos2/(Z_var2 + wgt2*median(Z_var2) + var_const2)^0.5
		}
	}

	A.max <- apply(abs(A), 2, max)
	B.max <- apply(abs(B), 2, max)
	tmp <- matrix(0,2, length(A.max))
	tmp[1,] <- A.max
	tmp[2,] <- B.max

	return(list(result = apply(tmp, 2, max), up.reg = A.max, down.reg = B.max, prof.up.reg = A, prof.down.reg = B))
}

###############################
# Now comes the sub-functions #
#                             #
###############################

GSZ_calc_prior <- function(var_table, mean_table_sq, num_genes, class_size){
	hyge.stat <- hyge.mean.var(num_genes, class_size)
	prob_sum <- count.prob.sum(num_genes, hyge.stat$mean, hyge.stat$var)
	Z_var1 <- 4*(var_table*prob_sum    + mean_table_sq*hyge.stat$var)
}

hyge.mean.var <- function(M,K){
	N <- c(1:M)
	out1 <- N*K/M
	out2 <- out1*((M-K)/M)*((M-N)/(M-1))
	return(list(mean= out1, var= out2))
}

count.prob.sum <- function(N, hyge.mean, hyge.var){
	N_tab <- c( 2:N )
	out1 <- hyge.mean[1]
	out2 <- (N_tab*hyge.mean[2:N] - hyge.mean[2:N]^2 - hyge.var[2:N])/(N_tab - 1) 
	return(c(out1, out2))
}






































