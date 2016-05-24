phylo.beta.pair<-function (x, tree, index.family = "sorensen")
{
    index.family <- match.arg(index.family, c("jaccard", "sorensen"))
   
    pbc<-x
    if (!inherits(x, "phylo.betapart")) {
    pbc <- phylo.betapart.core(x,tree)
    } # end of computing core results
	
	switch(index.family, sorensen = {
        phylo.beta.sim <- pbc$min.not.shared/(pbc$min.not.shared + pbc$shared)

        phylo.beta.sne <- ((pbc$max.not.shared - pbc$min.not.shared)/((2 * pbc$shared) + pbc$sum.not.shared)) * (pbc$shared/(pbc$min.not.shared + pbc$shared))

        phylo.beta.sor <- pbc$sum.not.shared/(2 * pbc$shared + pbc$sum.not.shared)

        phylo.pairwise <- list(phylo.beta.sim = phylo.beta.sim, phylo.beta.sne = phylo.beta.sne, phylo.beta.sor = phylo.beta.sor)
    								},

    					 jaccard = {
        phylo.beta.jtu <- (2 * pbc$min.not.shared)/((2 * pbc$min.not.shared) + pbc$shared)

        phylo.beta.jne <- ((pbc$max.not.shared - pbc$min.not.shared)/(pbc$shared + pbc$sum.not.shared)) * (pbc$shared/((2 * pbc$min.not.shared) + pbc$shared))

        phylo.beta.jac <- pbc$sum.not.shared/(pbc$shared + pbc$sum.not.shared)

        phylo.pairwise <- list(phylo.beta.jtu = phylo.beta.jtu, phylo.beta.jne = phylo.beta.jne, phylo.beta.jac = phylo.beta.jac)
    								}

    ) # end of switch

    return(phylo.pairwise)

} # end of function