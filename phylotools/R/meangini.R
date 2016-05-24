#### Function meangini as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

meangini <-
function(tree, subtree, times = 1000, plot = FALSE)
{
    tree$node.label <- NULL
    rH <- rh(tree, times = times)
    maxh <- max(as.hclust(multi2di(tree))$height)
    ineq <- c()
    for(i in 1:times){
        ineq[i] <- inequality(tree = tree, subtree = subtree, rH[i])
    }
    if(plot){
        par(mfrow = c(2,1))
        plot(tree)
        for(j in 1:times){
            axis(1)
            abline(v = (maxh)- rH[j], col = 2)
        }
        hist(ineq)
    }
    res <- mean(ineq)
    res <- list(meangini = res, pergini = ineq)
    return(res)
}

