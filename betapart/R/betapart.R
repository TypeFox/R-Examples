## multi site similarity values

beta.multi <- function(x, index.family="sorensen"){

	# test for a valid index
	index.family <- match.arg(index.family, c('jaccard','sorensen'))
	
	# test for pre-existing betapart objects
	if (! inherits(x, "betapart")){	
		x <- betapart.core(x)
	}

	maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
    minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])
	
	# run the analysis given the index
	switch(index.family,
		sorensen = {
            beta.sim <- minbibj / (minbibj + x$a)
            beta.sne <- (x$a / (minbibj + x$a)) * ((maxbibj - minbibj) / ((2 * x$a) + maxbibj + minbibj))
            beta.sor <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * x$a))

           	multi <- list(beta.SIM=beta.sim, beta.SNE=beta.sne,beta.SOR=beta.sor)},
		jaccard = {
            beta.jtu <- (2*minbibj) / ((2*minbibj) + x$a)
            beta.jne <- (x$a / ((2*minbibj) + x$a)) * ((maxbibj - minbibj) / ((x$a) + maxbibj + minbibj))
            beta.jac <- (minbibj + maxbibj) / (minbibj + maxbibj + x$a)

           	multi <- list(beta.JTU=beta.jtu, beta.JNE=beta.jne, beta.JAC=beta.jac)})

	return(multi)

}

## pairwise site similarity values

beta.pair <- function(x, index.family="sorensen"){
	
	# test for a valid index
	index.family <- match.arg(index.family, c('jaccard','sorensen'))
	
	# test for pre-existing betapart objects
	if (! inherits(x, "betapart")){	
		x <- betapart.core(x)
	}
	
	# run the analysis given the index
	switch(index.family,
		sorensen = {
			beta.sim <- x$min.not.shared / (x$min.not.shared + x$shared)
			beta.sne <- ((x$max.not.shared - x$min.not.shared) / ((2 * x$shared) + x$sum.not.shared)) * (x$shared / (x$min.not.shared + x$shared))
            beta.sor <- x$sum.not.shared / (2 * x$shared + x$sum.not.shared)
                
            pairwise <- list(beta.sim=as.dist(beta.sim), beta.sne=as.dist(beta.sne),beta.sor=as.dist(beta.sor))},
		jaccard = {
			beta.jtu <- (2*x$min.not.shared) / ((2*x$min.not.shared) + x$shared)
	        beta.jne <- ((x$max.not.shared - x$min.not.shared) / (x$shared + x$sum.not.shared)) * (x$shared / ((2*x$min.not.shared) + x$shared))
	        beta.jac <- x$sum.not.shared / (x$shared + x$sum.not.shared)

	        pairwise <- list(beta.jtu=as.dist(beta.jtu), beta.jne=as.dist(beta.jne),beta.jac=as.dist(beta.jac))})

	return(pairwise)
}