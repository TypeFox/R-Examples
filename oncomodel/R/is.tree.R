`is.tree` <-
function (tree) 
{
    if (!is.matrix(tree)) {
        stop("'tree' should be in 'matrix' format.")
        return(FALSE)
    }
    if (nrow(tree) != 2) {
        stop("'tree' should have two rows.")
        return(FALSE)
    }
    parents = match(tree[1, ], tree[2, ])
    if (any(tree[1, ] == tree[2, ])) {
        stop("tree has degenerate edges")
        return(FALSE)
    }
    if (length(which(parents < (1:length(parents)))) > 0) {
        cat("warning: tree has wrong order of edges\n")
        return(FALSE)
    }
    if (any(duplicated(tree[2, ]))) {
        cat("tree has non-unique parents\n")
        return(FALSE)
    }
    if (length(unique(tree[1, which(is.na(parents))])) > 1) {
        cat("tree has multiple roots\n")
        return(FALSE)
    }
    v <- sort(tree[2, !(1:ncol(tree) %in% parents)])
    if (!identical(v, as.numeric(1:length(v)))) {
        cat("leaves are not denoted by the first integers\n")
        return(FALSE)
    }
    return(TRUE)
}

