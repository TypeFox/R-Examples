SW <- function(n, nvertex, ncolor, edges, weights=1, beta){
    
    if (ncol(edges) != 2)
        stop("'edges' with two columns have to be provided.")
    nedge <- nrow(edges)
    if (nedge < length(weights))
        stop("The number of 'edges' is less than the number of 'weights'.")
    if (length(weights) < nedge){
            weights <- rep(weights, length=nrow(edges))
    }

    bondProbs <- 1 - exp(weights * (-beta))
    
    oneIteration <- sample(x=1:ncolor, nvertex, replace=TRUE) - 1
    oneIteration <- structure(as.integer(oneIteration), dim = dim(oneIteration))
    
    edges <- edges - 1
    edges <- structure(as.integer(edges), dim = dim(edges))

    colors <- .Call("sw", bondProbs, oneIteration, edges, nedge, as.integer(n), as.integer(nvertex), as.integer(ncolor))

    colors + 1
}

