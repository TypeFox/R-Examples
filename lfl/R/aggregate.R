aggregate <- function(conseq, degrees, partition) {
    if (!is.vector(conseq) || !is.character(conseq)) {
        stop("'conseq' must be a character vector")
    }
    if (!is.vector(degrees) || !is.numeric(degrees)) {
        stop("'degrees' must be a numeric vector")
    }
    if (length(conseq) != length(degrees)) {
        stop("The length of 'conseq' and 'degrees' must be the same")
    }
    if (!is.matrix(partition) || !is.numeric(partition)) {
        stop("'partition' must be a numeric matrix")
    }
    if (nrow(partition) <= 0 || ncol(partition) <= 0) {
        stop("'partition' must not be empty matrix")
    }
    if (is.null(colnames(partition)) || !is.character(colnames(partition))) {
        stop("'partition' must have colnames")
    }

    uniqConseq <- unique(conseq)
    if (length(intersect(uniqConseq, colnames(partition))) != length(uniqConseq)) {
        stop("Not all consequents are present in 'partition'")
    }

    res <- rep(1, nrow(partition))
    for (i in seq_along(conseq)) {
        part <- partition[, conseq[i]] + 1 - degrees[i]
        res <- pmin(res, part)
    }
    return(res)
}
