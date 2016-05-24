# get node names in pathways
pathway.nodes = function(pathway) {
    if(class(pathway) != "igraph") {
        stop("Wrong class for pathway.\n")
    }
    
    n = vcount(pathway)
    name = get.vertex.attribute(pathway, "name")
    
    if(length(name)) {
        return(name)
    }
    else {
        return(1:n)
    }
}
