#### Function inequality as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


inequality <- function(tree, subtree, h, detail = FALSE)
{
    a <- tree$tip.label
    b <- subtree$tip.label
    to.drop <- a[!a%in%b]
    tree$node.label <- NULL
    hclust.total <- hcreorder(as.hclust(multi2di(tree)))
    groups <- cutree(hclust.total, h = h)
    

    ngroups <- unique(groups)
    ndel <- c(); gn <- c()
    for(i in 1:length(ngroups)){
        groupi <- tree$tip.label[groups == i]
        gn[i] <- length(groupi)
        ## gn members for each group
        ndel[i] <- length(which(to.drop%in%groupi))
    }
    
    ratio <- ndel/gn
    datf <- data.frame(gn, ndel, ratio)
    giniindex <- gini(ratio)
    res = giniindex
    if(detail){
       res <- list(ratio, giniindex)
    }
    return(res)
}
