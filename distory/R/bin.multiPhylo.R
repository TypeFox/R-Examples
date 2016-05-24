bin.multiPhylo <- function(treelist)
{
    edge.dists <- dist.multiPhylo(treelist, method="edgeset")
    edge.dists <- as.matrix(edge.dists)

    bin.id = 1
    binning = rep(NA,length(treelist))
    while(any(is.na(binning)))
    {
         # find the next unbinned tree & find all matching trees
         next.tr = min(which(is.na(binning)))
         binning[edge.dists[,next.tr] == 0] = bin.id
         bin.id = bin.id + 1
    }

    binning
}
