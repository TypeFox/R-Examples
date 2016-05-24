#### Function resid.tree as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


resid.tree <-
function(tree, deltree){
    drop.tip(tree, deltree$retain)
}

