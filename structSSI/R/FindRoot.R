FindRoot <- function(tree.el) {
    tree.el.tmp <- tree.el
    root <- unique(tree.el.tmp[which(!(tree.el.tmp[, 1] %in% tree.el.tmp[, 2])), 1]) # no edge leads to root
    return(root)
}
