hamiltonian <- function(edges, start = 1, cycle = TRUE) {
    G <- graph_vect2list(edges)
    N <- length(G)
    X <- rep(0, N); X[1] <- start
    L <- vector("list", N); L[[1]] <- c()
    k <- 2; L[[2]] <- G[[start]]
    while (TRUE) {
        if (length(L[[k]]) > 0) {
            X[k] <- L[[k]][1]
            L[[k]] <- setdiff(L[[k]], X[k])
            if (k < N) {
                k <- k + 1
                L[[k]] <- setdiff(G[[X[k-1]]], X[1:(k-1)])
            } else {  # k == N
                if (!cycle || start %in% G[[X[k]]]) {
                    return(X)
                }
            }
        } else {
            if (k == 1) {
                return(NULL)
            } else {
                k <- k - 1
            }
        }
    }
}


## Graph as vector --> graph as list
graph_vect2list <- function(edges) {
    #   edges <- c(a1, a2, ...) for edges a1 -> a2, ...
    n <- length(edges)
    if (n %% 2 != 0 || any(edges <= 0) || any(floor(edges) != ceiling(edges)))
        stop("Argument 'edges' is not a correct edge list.")
    N <- max(edges)  # no. of vertices
    G <- vector(mode = "list", length = N)
    for (i in seq(1, length(edges), by = 2)) {
        v1 <- edges[i]; v2 <- edges[i+1]
        G[[v1]] <- c(G[[v1]], v2)
        G[[v2]] <- c(G[[v2]], v1)
    }
    for (i in 1:N) G[[i]] <- unique(G[[i]])
    return(G)
}
