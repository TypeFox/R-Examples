phylo.diff <- function(x, y, ...)
{
    uniqT1 <- distinct.edges(x, y)
    uniqT2 <- distinct.edges(y, x)
    
    treeA.cs <- rep("black", dim(x$edge)[1]) 
    treeA.cs[uniqT1] = "red"

    treeB.cs <- rep("black", dim(y$edge)[1]) 
    treeB.cs[uniqT2] = "red"

    par(mfrow=c(1,2))
    plot(x, edge.color=treeA.cs, ...)
    plot(y, edge.color=treeB.cs, ...)

    invisible()
}

distinct.edges <- function(x, y) # all edges in x not in y
{
    bp1 <- partition.leaves(x)
    bp1 <- lapply(bp1, sort)
    bp2 <- partition.leaves(y)
    bp2 <- lapply(bp2, sort)

    p = c()

    for(i in 1:length(bp1))
    {
        if(!(list(bp1[[i]]) %in% bp2))
        {
            p <- append(p, i)
        }
    }

    p
}

edge.from.split <- function(x, split)
{
    splits <- partition.leaves(x)
    splits <- lapply(splits, sort)
    split <- sort(split)
    match(list(split), splits)
}

get.bipartition <- function(x, e) # of leaves from edge
{
    acc = vector()
    inc = FALSE;
    for(i in 1:dim(x$edge)[1])
    {
        if(x$edge[i,1] == e)
        {
            acc <- append(acc, get.bipartition(x, x$edge[i,2]));
            inc = TRUE;
        }
    }

    if(!inc)
    {
        acc = x$tip.label[e];
    }

    acc
}

partition.leaves <- function(x) # get all bipartitoins
{
    r = list()
    for(i in 1:dim(x$edge)[1])
    {
        t <- get.bipartition(x, x$edge[i,2]);
        r <- append(r, list(t));
    }

    r
}

