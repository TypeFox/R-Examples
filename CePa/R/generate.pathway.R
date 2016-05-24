# construct the igraph object pathway
generate.pathway = function(el) {
    
    if(dim(el)[2] != 2) {
        stop("Second dimension of edgelist should be 2.\n")
    }
    
    # remove duplicat connections between two nodes
    el.vector = apply(el, 1, paste, collapse="--")
    el.vector = unique(el.vector)
    
    el = t(as.matrix(as.data.frame(strsplit(el.vector, "--"))))
    
    g = graph.edgelist(el, directed = TRUE)
    
    return(g)
}
