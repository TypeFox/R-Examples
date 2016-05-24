mult.choice <-
function (data, correct) {
    X <- data.matrix(data)
    p <- ncol(X)
    if (length(correct) != p)
        stop("'data' and 'correct' have incompatible dimensions.")
    ind <- vector("logical", p)
    for (i in 1:p)
        ind[i] <- correct[i] %in% unique(X[, i])
    if (!all(ind))
        stop("correct answer(s) for item(s) ", paste(which(!ind), collapse = ", "), " cannot be find in 'data'.")
    t(t(X) == correct) + 0.0
}
