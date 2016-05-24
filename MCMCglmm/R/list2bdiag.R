"list2bdiag"<-
function (x)
{
    d <- matrix(unlist(lapply(x, dim)), ncol = 2, byrow = TRUE)
    if (nrow(d) != length(x))
        stop("all arguments must be matrices")
    dsum <- apply(d, 2, sum)
    z <- array(0, dsum)
    dimnames(z) <- list(character(dsum[1]), character(dsum[2]))
    coord <- c(0, 0)
    for (i in 1:length(x)) {
        coord1 <- coord[1] + 1:d[i, 1]
        coord2 <- coord[2] + 1:d[i, 2]
        z[coord1, coord2] <- x[[i]]
        rownames(z)[coord1] <- rownames(x[[i]], do.NULL = FALSE,
            prefix = "")
        colnames(z)[coord2] <- colnames(x[[i]], do.NULL = FALSE,
            prefix = "")
        coord <- coord + d[i, ]
    }
    z
}
