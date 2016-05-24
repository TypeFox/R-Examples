load_gct <- function(file = NULL)
{
	gct_info <- readLines(file)
	col_names <- noquote(unlist(strsplit(gct_info[3], "\t")))
	col_names <- col_names[-c(1:2)]
	gct_info <- gct_info[-c(1:3)]
	n_gene <- length(gct_info)
	n_sample <- length(col_names)
	row_names <- vector(length = n_gene, mode = "character")
	gene_anno <- vector(length = n_gene, mode = "character")
	expmat <- matrix(0, nrow = n_gene, ncol = n_sample)
	for(i in 1:n_gene){
		tmp <- noquote(unlist(strsplit(gct_info[i], "\t")))
		row_names[i] <- tmp[1]
		gene_anno[i] <- tmp[2]
		tmp <- tmp[-c(1:2)]
		mode(tmp) <- "numeric"
		expmat[i, ] <- tmp
	}
	rownames(expmat) <- row_names
	colnames(expmat) <- col_names
	return(expmat)
}
