exogenize <- function(mat, cmx, items = seq_len(ncol(mat)), endnode, crossitem = NULL) {
    stopifnot(is.matrix(mat),
              length(unique(as.vector(mat))) == 2, # must be a binary matrix
              is.integer(items),
              min(items) >= 1,
              max(items) <= ncol(mat),
              length(unique(items)) == (li <- length(items)), # no repetitions
              length(endnode <- factor(endnode)) == li,
              is.matrix(cmx),
              (nr <- nrow(cmx)) == length(levels(endnode)))
    nc <- ncol(cmx)
    m1 <- droplevels(subset(tolong(mat), as.integer(item) %in% items))
    intItem <- as.integer(m1$item)
    m1$endnode <- endnode[intItem]
    if (!is.null(crossitem)) {
        stopifnot(length(crossitem <- factor(crossitem)) == li)
        m1$crossitem <- crossitem[intItem]
    }
    colnames(cmx) <- sprintf(paste("exo%0", nchar(nc), "d", sep=''), seq_len(nc))
    cbind(m1, cmx[as.integer(m1$endnode),])
}
