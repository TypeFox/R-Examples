### Add lines from vegan::spantree object to an ordirgl plot
`orglspantree` <-
    function(object, spantree, display = "sites", choices = 1:3,
             col = "black", ...)
{
    if (!inherits(spantree, "spantree"))
        stop("'spantree' must be a result of vegan::spantree() function")
    x <- scores(object, display = display, choices = choices, ...)
    ## get kids
    k <- spantree$kid
    ## change colors to rgb in 0..1 and recycle if needed
    col <- rep(col, length = nrow(x))
    if (is.factor(col))
        col <- as.numeric(col)
    col <- col2rgb(col)/255
    ## lines
    for (i in 1:length(k)) {
        if (is.na(k[i])) # skip NA links: disconnected spantree
            next
        one <- x[i+1,]
        two <- x[k[i],]
        lcol <- rgb(t(col[, i+1] + col[,k[i]])/2)
        rgl.lines(rbind(one, two), col = lcol, ...)
    }
}

### Add lines from an hclust object to an ordirgl plot

`orglcluster` <-
    function(object, cluster, prune = 0, display = "sites",
             choices = 1:3, col = "black", ...)
{
    if (!inherits(cluster, "hclust"))
        cluster <- as.hclust(cluster)
    x <- scores(object, display = "sites", choices = choices, ...)
    ## recycle colours if needed
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = nrow(x))
    ## (Ab)use vegan:::reorder.hclust to get the coordinates and
    ## colours of internal nodes
    node <- apply(x, 2, function(val) reorder(cluster, val)$value)
    nodecol <- apply(col2rgb(col)/255, 1,
                     function(val) reorder(cluster, val)$value)
    nodecol <- rgb(nodecol)
    ## go through merge matrix
    merge <- cluster$merge
    for(i in seq_len(nrow(merge) - prune)) {
        if(merge[i,1] < 0)
            one <- x[-merge[i,1],]
        else
            one <- node[merge[i,1],]
        if (merge[i,2] < 0)
            two <- x[-merge[i,2],]
        else
            two <- node[merge[i,2],]
        rgl.lines(rbind(one, two), col = nodecol[i], ...)
    }
}
