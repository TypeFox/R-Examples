#### Function del.tree.tip as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

del.tree.tip <-
function(tree, n)
{
    N <- length(tree$tip.label)
    tip.nam <- sample(1:N, size = n)
    to.drop <- (tree$tip.label)[tip.nam]
    retain  <- (tree$tip.label)[-tip.nam]
    #to.drop <- c("Struthioniformes", "Tinamiformes", "Craciformes")
    subtree <- drop.tip(tree, to.drop)
    list(subtree = subtree, to.drop = to.drop, retain = retain)
}

