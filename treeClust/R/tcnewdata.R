tcnewdata <- function (obj, d.num = 1, tbl, mat, trees) 
{
#
# Compute the "newdata" data frame from a treeClust item. If there's
# a newdata element in there, computed with the proper d.num, return it.
#
if (!missing (obj) && any(names (obj) == "newdata") && obj$d.num == d.num)
    return (obj$newdata)
if (!missing (obj) && any (names (obj) == "tbl"))
    tbl <- obj$tbl
else
    if (missing (tbl)) 
        if (d.num == 2 || d.num == 4)
            stop ("'Tbl' missing but required")

if (!missing(obj) && any(names(obj) == "mat")) 
    mat <- obj$mat
else if (missing(mat)) 
    stop("'Mat' misisng but required")

#
# We need to know the number of leaves in each tree (from which we subtract
# 1 to get the number of columns) and the sample size.
# Then we can set up newdata.
#
leaf.counts <- sapply(mat, function(x) length(unique(x)))
n <- nrow (mat)
start <- c(1, 1 + cumsum(leaf.counts[-length(leaf.counts)]))
end <- cumsum(leaf.counts)
newdata <- matrix(0, n, sum(leaf.counts))

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

    for (i in 1:length(leaf.counts)) {
        mod <- model.matrix(~factor(mat[, i]) -  1)
        if (d.num == 2) 
          newdata[, start[i]:end[i]] <- mod * tbl[i, 
            "DevRat"]/max(tbl[, "DevRat"])
        else newdata[, start[i]:end[i]] <- mod
    }
    return (newdata)
}
#
# For d3 or d4 we need the trees. Also we only produce (# leaves - 1)
# columns for each tree.
#
col.counts <- sapply(mat, function(x) length(unique(x)) - 1)
n <- nrow (mat)
start <- c(1, 1 + cumsum(col.counts[-length(col.counts)]))
end <- cumsum(col.counts)
newdata <- matrix(0, n, sum(col.counts))

if (!missing (obj) && any (names (obj) == "trees")) 
    trees <- obj$trees
else
    if (missing (trees))
        stop ("For d3 or d4, we need the trees")
#
# Compute the matrix of pairwise leaf distances, then call
# cmdscale() on the result.
#
for (i in 1:length (trees)) {
    leaf.dists <- d3.dist (trees[[i]], return.pd=TRUE)
    newcols <- cmdscale (leaf.dists, ncol(leaf.dists) - 1)
#
# The row.names of "newcols" are the leaf numbers, not the numbers
# found in the "where" element of the tree. So we convert...
#
    w <- trees[[i]]$where
    r <- row.names (trees[[i]]$frame)[w]
#
# ...and then extract the correct rows of newcols.
#
    if (d.num == 4) 
            newdata[, start[i]:end[i]] <- newcols[r,] * tbl[i, 
            "DevRat"]/max(tbl[, "DevRat"])
        else  {

            newdata[, start[i]:end[i]] <- newcols[r,]
        }
}
return (newdata)
}

