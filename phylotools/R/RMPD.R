#### Function RMPD as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


RMPD <-
function(subtree, tree){
    mean(as.dist(cophenetic(subtree)))/mean(as.dist(cophenetic(tree)))
}

