phylo.beta.multi<-function (x, tree, index.family = "sorensen")
{
    index.family <- match.arg(index.family, c("jaccard", "sorensen"))
    pbc<-x
    if (!inherits(x, "phylo.betapart")) {
    pbc <- phylo.betapart.core(x,tree)
    } # end of computing core results

    switch(index.family, sorensen = {
        phylo.beta.SIM <- sum(pbc$min.not.shared)/(pbc$sumSi - pbc$St + sum(pbc$min.not.shared))

        phylo.beta.SNE <- ((sum(pbc$max.not.shared) - sum(pbc$min.not.shared))/(2 * (pbc$sumSi - pbc$St) + sum(pbc$min.not.shared) + sum(pbc$max.not.shared))) * ((pbc$sumSi - pbc$St)/(pbc$sumSi - pbc$St + sum(pbc$min.not.shared)))

        phylo.beta.SOR <- (sum(pbc$min.not.shared) + sum(pbc$max.not.shared))/(2 * (pbc$sumSi - pbc$St) + sum(pbc$min.not.shared) + sum(pbc$max.not.shared))

        phylo.multi <- list(phylo.beta.SIM = phylo.beta.SIM, phylo.beta.SNE = phylo.beta.SNE, phylo.beta.SOR = phylo.beta.SOR)
    								},

    					 jaccard = {
        phylo.beta.JTU <- (2 * sum(pbc$min.not.shared))/((2 * sum(pbc$min.not.shared)) + pbc$sumSi - pbc$St)

        phylo.beta.JNE <- ((sum(pbc$max.not.shared) - sum(pbc$min.not.shared))/(pbc$sumSi - pbc$St + sum(pbc$max.not.shared) + sum(pbc$min.not.shared))) * ((pbc$sumSi - pbc$St)/(2 * sum(pbc$min.not.shared) + pbc$sumSi - pbc$St))

        phylo.beta.JAC <- (sum(pbc$min.not.shared) + sum(pbc$max.not.shared))/(pbc$sumSi - pbc$St + sum(pbc$min.not.shared) + sum(pbc$max.not.shared))

        phylo.multi <- list(phylo.beta.JTU = phylo.beta.JTU, phylo.beta.JNE = phylo.beta.JNE, phylo.beta.JAC = phylo.beta.JAC)
    								}

    ) # end of switch

    return(phylo.multi)

} # end of function