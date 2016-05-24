##Script generated in:
# 2011
# 10:39:53 AM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

correspondence <- function(index, sign, skein.sign) {
	index <- index + 1
	p <- index * sign
	set <- c(1, 2, -1, -2)
	p <- which(p == set)
	correspondence <- switch(p,
			paste(-skein.sign, "*l**-2", sep = ""), ## 1
			paste(-skein.sign, "*l**-1*m", sep = ""), ## 2
			paste(-skein.sign, "*l**2", sep = ""), ## -1
			paste("l*-m", sep = ""), ## -2
	)
	return(correspondence)
} 


contribute <- function(ancestor, signs, skein.sign) {
	n <- length(ancestor)
	single.contr <- rep(Inf, n)
	for(i in 1 : n) {
		single.contr[i] <- correspondence(ancestor[i], signs[i], skein.sign)
	}
	global.contr <- paste("Mul(", paste(single.contr, collapse = ",") ,")")
	return(global.contr)
}

HOMFLYpolynomial <- function(leaves, tree, skein.sign  = -1) {
	
	sympy(Var('l'))
	sympy(Var('m'))
	
	if(identical(leaves, list())) {
		nends <- length(tree[[1]]$ends)
		toeval <- paste('(m**-1 *', '(', -skein.sign, '*l-l**-1))**', nends)
		polynomial <- sympy(toeval)
	}
	else {
		n <- length(leaves)
		tree.indices <- lapply(leaves, bin2N)
		nendsv <- lapply(tree.indices, function(x) length(tree[[x]]$ends))		
		ancestors <- lapply(leaves, function(x) x[-1])
		signs <- lapply(tree.indices, function(x) tree[[x]]$signs)
		ancestors.contribute <- rep(Inf, n)
		components.contribute <- rep(Inf, n)
		for(i in 1 : n) {
			components.contribute[i] <- 
					paste('(m**-1 *', '(', -skein.sign, '*l-l**-1))**', nendsv[[i]])
			ancestors.contribute[i] <- contribute(ancestors[[i]], signs[[i]], skein.sign)
		}
		toeval <- paste('Add(', 
				paste(components.contribute, '*', ancestors.contribute, collapse = '+'), ').expand()')
		#print(toeval)
		polynomial <- sympy(toeval)
		polynomial <- parseToR(polynomial)
	}
	return(polynomial)
}

