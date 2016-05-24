tcdist <- function (obj, d.num = 1, tbl, mat, trees, verbose = 0)
{
#
# Extract distances from an existing treeClust object. If there's
# a dist element in there, computed with the proper d.num, return it.
#
if (!missing (obj) && any(names (obj) == "dists") && obj$d.num == d.num)
    return (obj$dists)
if (!missing (obj) && any (names (obj) == "tbl"))
    tbl <- obj$tbl
else
    if (missing (tbl)) 
        if (d.num == 2 || d.num == 4)
            stop ("'Tbl' missing but required")
#
# For d1, use daisy() on the "mat" object; for d2, do the same
# but with weights.
#
if (d.num == 1 || d.num == 2) {
    if (!missing (obj) && any (names (obj) == "mat"))
        mat <- obj$mat
    else
        if (missing (mat))
            stop ("For d1 or d2, this function requires the 'mat' element")
    if (d.num == 1) tree.wts <- rep (1, ncol(mat))
    if (d.num == 2) tree.wts <- tbl[,"DevRat"] / max (tbl[,"DevRat"])
    return (daisy (mat, metric = "gower", weights = tree.wts))
}
#
# For d3 or d4 we need the trees.
#
if (!missing (obj) && any (names (obj) == "trees")) 
    trees <- obj$trees
else
    if (missing (trees))
        stop ("For d3 or d4, we need the trees")
n <- length (trees[[1]]$where)
dists <- numeric (n * (n - 1) / 2)
if (d.num == 3) tree.wts <- rep (1, length(trees))
if (d.num == 4) tree.wts <- tbl[,"DevRat"] / max (tbl[,"DevRat"])

for (i in 1:length (trees)) {
    if (verbose > 0)
        cat ("Tree ", i, ", has wt ", tree.wts[i], "\n")
    dists <- dists + tree.wts[i] * d3.dist (trees[[i]])
}
class (dists) <- "dist"
attr (dists, "Size") <- n
attr (dists, "Diag") <- FALSE
attr (dists, "Upper") <- FALSE
attr (dists, "method") <- "manhattan"
return (dists)
}

