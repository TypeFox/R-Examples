ListTreePval <- function(tree) {
    # Given an igraph tree, returns format for d3Network
    # Assumes tree has V(tree)$names equal to names we want
    # Also assumes tree has p-values associated with it
    if(length(V(tree)) == 1) {
        return(list(name = V(tree)$name, pval = V(tree)$pval))
    } else {
        parent <- V(tree)[1]$name
        cur.pval <- V(tree)[1]$pval
        edgelist <- get.edgelist(tree, names = FALSE)
        root.position.in.edgelist <- which(edgelist[,1] == 1)
        children <- edgelist[root.position.in.edgelist, 2]
        children.list <- list()
        for(i in 1:length(children)) {

            subcomp.indices <- subcomponent(tree, children[i], 'out')
            subgraph <- induced.subgraph(graph = tree,
                                         vids = subcomp.indices)
            children.list[[i]] <- ListTreePval(induced.subgraph(graph = tree,
                                                            vids = subcomp.indices))
        }
        result <- list(name = parent, pval = cur.pval, children = children.list)
    }
    return (result)
}
