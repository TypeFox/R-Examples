`permute` <- function(i, n, control) {
    complete <- getComplete(control)
    ap <- getAllperms(control)
    perm <- if (complete && !is.null(ap)) {
        ap[i, ]                 # select ith permutation
    } else {
        if (complete) {
            warning("'$all.perms' is NULL, yet '$complete = TRUE'.\nReturning a random permutation.")
        }
        shuffle(n, control)
    }
    perm
}
