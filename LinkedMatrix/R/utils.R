chunkRanges <- function(a, n, i = NULL) {
    if (n > a) {
        stop(paste("Cannot split", a, "into", n, "chunks. Reduce the number of chunks."))
    }
    k <- as.integer(a / n)
    r <- as.integer(a %% n)
    range <- function(i, k, r) {
        c((i - 1) * k + min(i - 1, r) + 1, i * k + min(i, r))
    }
    if (!is.null(i)) {
        range(i, k, r)
    } else {
        sapply(seq_len(n), range, k, r)
    }
}
