treeSims<-function (node.height,perms=100, tip.label = NULL, br = "coalescent", ...) 
{
#this is slightly modified version of function 'rcoal' from APE	

	node.height<-sort(node.height)
	n <- length(node.height)+1
    nbr <- 2 * n - 2
    edge <- matrix(NA, nbr, 2)
    x <- node.height
    res<-list()
    for (j in 1:perms){
    

        edge.length <- numeric(nbr)
        h <- numeric(2 * n - 1)
        #node.height <- cumsum(x)
        pool <- 1:n
        nextnode <- 2L * n - 1L
        for (i in 1:(n - 1)) {
            y <- sample(pool, size = 2)
            ind <- (i - 1) * 2 + 1:2
            edge[ind, 2] <- y
            edge[ind, 1] <- nextnode
            edge.length[ind] <- node.height[i] - h[y]
            h[nextnode] <- node.height[i]
            pool <- c(pool[!pool %in% y], nextnode)
            nextnode <- nextnode - 1L
        }
#    }
    phy <- list(edge = edge, edge.length = edge.length)
    if (is.null(tip.label)) 
        tip.label <- paste("t", 1:n, sep = "")
    phy$tip.label <- sample(tip.label)
    phy$Nnode <- n - 1L
    class(phy) <- "phylo"
    phy <- reorder(phy)
    phy$edge[phy$edge[, 2] <= n, 2] <- 1:n
   res[[j]]<-phy
   }
   return(res)
}
