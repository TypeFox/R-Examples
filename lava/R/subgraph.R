subgraph <- function(g,from,to,Tree=new("graphNEL",node=c(to,from),edgemode="directed"),...) {
    adjnodes <- graph::adj(g,from)[[1]]
    newnodes <- !(adjnodes %in% graph::nodes(Tree))
    if (length(adjnodes)==0)
        return(Tree)
    for (v in adjnodes) {
        if (v==to) {
            Tree <- graph::addEdge(from, v, Tree)
        }
        re1 <- graph::acc(g,v)[[1]] ## Reachable nodes from v
        if ((to %in% names(re1)[re1>0])) {
            if (!(v %in% graph::nodes(Tree)))
                Tree <- graph::addNode(v,Tree)
            Tree <- graph::addEdge(from, v, Tree)
            Tree <- path(g,v,to,Tree)
        }
    }
    return(Tree)
}
