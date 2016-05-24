`dino.mst` <-
function (x, random.start = TRUE, random.search = TRUE) {
    X <- as.matrix(x)
    n <- dim(X)[1]
    N <- matrix(0, n, n)
    large.value <- max(X) + 1
    diag(X) <- large.value
    if (random.start == TRUE) 
        tree <- sample(1:n, 1)
    else tree <- 1
    ntree <- c(1:n)[-tree]
    while (length(tree) < n) {
        m <- min(X[ntree, tree])
        ind <- which(as.matrix(X[ntree, tree]) == m, TRUE)
        li <- length(ind[, 1])
        if (li > 1 & random.search == TRUE) no <- sample(1:li, 1)
        else no <- 1
        if (length(ntree) == 1) {
            nti <- ntree[ind[no, 2]]
            ti <- tree[ind[no, 1]]
        }
        else {
            nti <- ntree[ind[no, 1]]
            ti <- tree[ind[no, 2]]
        }
        N[nti, ti] <- 1
        N[ti, nti] <- 1
        cs <- colSums(N)
        tree <- which(cs > 0)
        ntree <- c(1:n)[-tree]
    }
    dimnames(N) <- dimnames(X)
    return(N)
}

