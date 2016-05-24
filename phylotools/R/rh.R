#### Function rh as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

rh <-
function(tree, times = 1)
{
    tree$node.label <- NULL
    h = (as.hclust(multi2di(tree))$height)
    hmax = max(h)
    hmin = min(h)
    res <- runif(hmin, hmax,n = times)
    return(res)
}

