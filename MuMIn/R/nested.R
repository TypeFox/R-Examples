`nested` <- 
function(x, indices = c("none", "numeric", "rownames"), rank = NULL) {

	indices <- match.arg(indices)

	if(!inherits(x, "model.selection")) 
		stop("'x' is not a \"model.selection\" object")
		
	vColIdx <- type2col(x, "varying")
	if(nVCols <- length(vColIdx)) {
		vtab <- x[, vColIdx, drop = FALSE]
		for(i in 1L:ncol(vtab)) vtab[, i] <- as.numeric(vtab[, i])
		vtab <- as.matrix(vtab)
	}	
		
	tab <- !is.na(x[, attr(x, "terms")])
	n <- nrow(tab)

	if(indices == "none") {
		if(is.null(rank)) {
			rank <- colnames(x)[which(colnames(x) == "delta")[1L] - 1L]
		} else if (!is.na(rank) && !rank %in% colnames(x))
			cry(, "column named \"%s\" does not exist in 'x'", rank)
		
		if(!is.na(rank) && any(diff(x[, rank]) < 0))
			cry(, "'x' is not ordered by \'%s\'", rank, warn = TRUE)
	
		
		is.nested <- function(x, inside) all(inside == x | inside)
		vmatch <- if(nVCols) 
			function(i, j) vtab[i, ] == vtab[j, ] else
			function(i, j) TRUE
			
		res <- logical(n)
		for(i in 2L:n)
			for(j in seq_len(i - 1L))
				if(vmatch(i, j) && is.nested(tab[j, ], tab[i, ])) {
					res[i] <- TRUE
					break;
				}
	} else {
	
		# 'alldf': same as apply(g, margin, all) but ~2x faster
		alldf <- function(g, margin = 1L) {
			g <- as.matrix(g)
			dg <- dim(g)
			mode(g) <- "integer"
			if(margin == 1L) .rowSums(g, dg[1L], dg[2L]) == dg[2L] else
				.colSums(g, dg[1L], dg[2L]) == dg[1L]
		}
	
		vmatch <- if(nVCols) 
			function(i) alldf(apply(vtab, 1L, "==", vtab[i, ]), 1L) else
			function(i) TRUE
		tab <- t(tab)
		idx <- seq.int(n)
		res <- vector(length = n, mode = "list")
		for(i in idx) {
			z <- vmatch(i) & alldf(tab == tab[, i] | tab[, i], 2L)
			z[i] <- FALSE
			res[[i]] <- which(z)
		}
		if(indices == "rownames") {
			res <- lapply(res, function(i, x) x[i], rownames(x))
		}
		names(res) <- rownames(x)
	}
	return(res)
}
