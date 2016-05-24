tclos <-
function (X, labels = NULL) 
{
    Y <- as.matrix(X)
    Y <- replace(Y, Y == "p", 1)
    Y <- replace(Y, Y != 1, 0)
    for (i in seq_len(ncol(Y))) {
        tmp <- outer(Y[, i], Y[i, ], pmin.int)
        Y <- pmax(Y, tmp)
    }
    Y <- replace(Y, Y == 1, "p")
    Y <- replace(Y, Y != "p", "o")
    n <- which(X == "n", arr.ind = TRUE)
    for (i in 1:nrow(n)) Y[n[, 1][i], n[, 2][i]] <- "n"
    rm(i, n)
    a <- which(X == "a", arr.ind = TRUE)
    for (i in 1:nrow(a)) Y[a[, 1][i], a[, 2][i]] <- "a"
    rm(i, a)
    q <- which(X == "q", arr.ind = TRUE)
    for (i in 1:nrow(q)) Y[q[, 1][i], q[, 2][i]] <- "q"
    rm(i, q)
    X <- data.frame(matrix(nrow = nrow(X), ncol = ncol(X)))
    for (i in 1:nrow(Y)) X[i, ] <- Y[i, ]
    rm(i)
    if (isTRUE(is.null(labels) == FALSE) == TRUE) 
        rownames(X) <- colnames(X) <- labels
    X
}
