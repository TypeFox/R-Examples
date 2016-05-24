#### Function phyloshuffle as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

phyloshuffle <- function(tree)
{
    labels <- tree$tip.label
    randmz.label <- sample(labels)
    tree$tip.label <- randmz.label
    return(tree)
}

