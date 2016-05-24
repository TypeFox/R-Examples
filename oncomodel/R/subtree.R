`subtree` <-
function (node, tree) 
{
    edges <- numeric(0)
    newedges <- which(tree[1, ] == node)
    while (length(newedges) > 0) {
        edges <- c(edges, newedges)
        newnodes <- tree[2, newedges]
        newedges <- which(tree[1, ] %in% newnodes)
    }
    return(sort(edges))
}

