convert.graph <- function(graph) {
    if (is.matrix(graph))
        t(graph)
    else if (is.data.frame(graph))
        t(data.matrix(graph))
    else if (inherits(graph, "graph") && requireNamespace("graph", quietly=TRUE)) {
        graph::edgeMatrix(graph)
    }
    else stop("unrecognized graph type")
}

name.orbits <- function(orbits) {
    orb.names = NULL
    for (i in 0:(ncol(orbits)-1)) {
        orb.names <- c(orb.names, paste("O", i, sep=""))
    }
    colnames(orbits) <- orb.names
    orbits
}

count4 <- function(graph) {
    edges <- convert.graph(graph) 
    result <- .C("count4",
        edges, dim(edges),
        orbits=matrix(0, nrow=max(edges), ncol=15))$orbits
    name.orbits(result)
}

count5 <- function(graph) {
    edges <- convert.graph(graph) 
    result <- .C("count5",
        edges, dim(edges),
        orbits=matrix(0, nrow=max(edges), ncol=73))$orbits
    name.orbits(result)
}

ecount4 <- function(graph) {
    edges <- convert.graph(graph) 
    result <- .C("ecount4",
        edges, dim(edges),
        orbits=matrix(0, nrow=ncol(edges), ncol=12))$orbits
    name.orbits(result)
}

ecount5 <- function(graph) {
    edges <- convert.graph(graph) 
    result <- .C("ecount5",
        edges, dim(edges),
        orbits=matrix(0, nrow=ncol(edges), ncol=68))$orbits
    name.orbits(result)
}
