`comp.freq` <-
function (x, tree, p) 
{
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("'x' should be in 'matrix' format.")
    if (!is.tree(tree)) 
        stop("'tree' is not in suitable matrix format.")
    if (!is.vector(p) | length(p) != ncol(tree)) 
        stop("Sizes of 'tree' and 'p' are not compatible.")
    singlefreq <- apply(x, 2, mean)
    pairfreq <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    for (j in 1:ncol(x)) {
        for (k in 1:ncol(x)) pairfreq[j, k] <- (x[, j] %*% x[, 
            k])/nrow(x)
    }
    singleprob <- rep(1, ncol(tree))
    for (j in tree[2, ]) {
        node <- j
        while (any(tree[2, ] == node)) {
            edge <- which(tree[2, ] == node)
            singleprob[j] <- singleprob[j] * p[edge]
            node <- tree[1, edge]
        }
    }
    if (interactive()) 
        x11(width = 8, height = 4)
    par(mfrow = c(1, 2))
    plot(singleprob[1:ncol(x)], singlefreq, main = "single probabilities vs. frequencies")
    pairprob <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    for (j in 1:ncol(x)) {
        for (k in 1:ncol(x)) pairprob[j, k] <- singleprob[j] * 
            singleprob[k]/singleprob[mrca(j, k, tree)]
    }
    plot(pairprob, pairfreq, main = "pair probabilities vs. frequencies")
}

