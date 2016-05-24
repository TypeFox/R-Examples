library(phylobase)
library(ape)

set.seed(1)
r1 <- rcoal(5)

## trace("phylo4d", browser, signature = "phylo")
## untrace("phylo4d", signature = "phylo")
tipdat <- data.frame(a=1:5,row.names=r1$tip.label)
p1 <- phylo4d(r1,tip.data=tipdat,node.data=data.frame(a=6:9), match.data=FALSE)
p2 <- prune(p1,1)
summary(p2)

## from picante
`phylo2phylog` <-
function(phy, ...) {
    newick2phylog(write.tree(phy, multi.line = FALSE),...)
}

plot.phylo(as(p2,"phylo"))
