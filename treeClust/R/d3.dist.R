d3.dist <- function (mytree, return.pd=FALSE) {
#
# d3 distance. Given a tree, produce a set of columns that describe
# the d3-style inter-point distances. Leaf numbers are the usual node 
# numbers; leaf "F.nums" are the row numbers within "frame", as recorded
# in the "where" element.
#
myframe <- mytree$frame
#
# Function to find the nearest ancestor to two leaves. "Paths" comes
# from make.leaf.paths.
#
leaf.ancestor <- function (L1, L2, paths)
{
    out <- numeric (length (L1))
    for (i in 1:length (L1)) {
        path1 <- paths[L1[i],]
        path2 <- paths[L2[i],]
        p1 <- path1[path1 != 0]
        p2 <- path2[path2 != 0]
        out[i] <- max (intersect (p1, p2))
    }
    return (out)
}
#
# Funtion to identify leaves that are children of row "num". If num
# is odd, that's all the leaves with larger numbers, lower down in the
# "myframe" than "num," until we encounter a smaller number. (E.g. if num
# is 11, we get all the leaves below us until we run into node 3, 6, or 7.)
# if num is even, we stop at num+1 -- e.g. if num is 10, we stop at nodes
# 3, 6, 7 or 11.
#
leaves.beneath.num <- function (num, myframe) {
    if (num == 1) return (row.names(myframe)[myframe[,"var"] == "<leaf>"])
    myrow <- which (row.names(myframe) == paste (num))
    if (length (myrow) == 0) return (0)
    upper.limit <- ifelse (num %% 2 == 0, num + 1, num)
    smaller.numbers <- as.numeric (row.names (myframe)) <= upper.limit
    if (!any (smaller.numbers[(myrow+1):nrow(myframe)])) {
        subset <- myframe[(myrow+1):nrow (myframe),]
    } else {
        higher.rows <- (1:nrow(myframe)) > myrow
        if (!any (smaller.numbers & higher.rows)) return (0)
        one.too.far <- min (which (smaller.numbers & higher.rows))
        subset <- myframe[myrow:(one.too.far-1),]
    }
    return (row.names (subset)[subset[,"var"] == "<leaf>"])
}
#
leaf.nums <- as.numeric (row.names (myframe[myframe[,"var"] == "<leaf>",]))
leaf.paths <- make.leaf.paths (max (leaf.nums))
leaf.F.nums <- which (myframe[,"var"] == "<leaf>")
#
# Figure out the immediate ancestor for each pair of leaves.
#
ancestors <- outer (leaf.nums, leaf.nums, leaf.ancestor, paths = leaf.paths)
dimnames (ancestors) <- list (leaf.nums, leaf.nums)
#
# Get the set of interior nodes. 
# The set of leaves descendant from an interior 
# node i is the set beneath it on the list, up to the first time we encounter 
# a node number smaller than (i+1), if i is even, or i, if i is odd, if 
# there is one. 
#
intnode.nums   <- as.numeric (row.names (myframe[myframe[,"var"] != "<leaf>",]))
intnode.F.nums <- which (myframe[,"var"] != "<leaf>")

#-------------------------
deviances <- rp.deviance (mytree)
intnode.sums <- numeric (length (intnode.nums))
for (i in 1:length (intnode.sums)) {
    lbn <- leaves.beneath.num (intnode.nums[i], myframe)
    intnode.sums[i] <- sum (deviances[is.element (row.names (myframe), lbn)])
}
names (intnode.sums) <- intnode.nums
#
# Intnode sums [i] gives the sum of the deviance in the leaves
# associated with pruning back to leaf i. The cost to us of doing
# that pruning is (dev. at this leaf - that sum).
#
intnode.dev <- deviances[intnode.F.nums]
root.dev <- deviances[1]
leaf.dev <- sum (deviances[leaf.F.nums])

dev.dist <- (intnode.dev - intnode.sums) / (root.dev - leaf.dev)
#
# Find the ancestor for each pair of leaves.
#
pairwise.distances <- ancestors
pairwise.distances[] <- dev.dist[as.character(ancestors)]
diag (pairwise.distances) <- 0
if (return.pd)
    return (pairwise.distances)
#
# The "where" item identifies each observation by its row number
# in the frame, not by the leaf number. So convert it to leaf number.
#
where <- as.numeric (row.names (mytree$frame)[mytree$where])
#
# Now, the leaves index the rows and columns of pairwise.distances.
# So convert from leaf number to row-of-pairwise-distance.
#
rowholder <- match (where, row.names (pairwise.distances))
rn <- length(where)
#
# We could expand.grid(), but that would create n^2 pairs, and
# we'd like to keep that number down. So we'll construct them ourselves.
# 
grid <- matrix (0, nrow = rn * (rn - 1)/2, ncol=2)
grid[,1] <- rep (1:(rn - 1), (rn - 1):1)
grid[,2] <- unlist (sapply (2:rn, function (x) seq (x, rn)))
dists <- pairwise.distances[cbind(rowholder[grid[, 1]], rowholder[grid[, 
        2]])]
class(dists) <- "dist"
attr(dists, "Size") <- length (mytree$where)
attr(dists, "Diag") <- FALSE
attr(dists, "Upper") <- FALSE
attr(dists, "method") <- "manhattan"
return(dists)
}


