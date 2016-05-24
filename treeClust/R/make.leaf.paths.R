make.leaf.paths <- function (up.to = 2047) 
{
#
# Create a matrix of leaf paths. This has one row for each integer, and
# the values in that row give the parents of that leaf number. For example,
# row 122 has 1, 3, 7, 15, 30, and 61. The largest number of columns
# possible is log-base-2 of up.to, rounded up to the next integer. Exception:
# if up.to is exactly a power of 2, add 1. For a tree with two leaves,
# just build and return the thing.
#
if (up.to == 3)
     return (matrix (c(1, 1, 1, 0, 2, 3), ncol=2))
num.of.cols <- ceiling (log (up.to, 2))
if (2^num.of.cols == up.to) num.of.cols <- num.of.cols + 1
mat <- matrix (0, nrow = up.to, ncol = num.of.cols)
row.names (mat) <- 1:up.to
#
# Do this in a loop. It's not efficient, but we only do this once.
# First of all, place the values for 1, 2, and 3.
#
mat[1,1] <- 1
mat[2,1:2] <- c(1, 2)
mat[3,1:2] <- c(1, 3)
#
# For each subsequent row, grab its parent and append the current index.
#
for (i in 4:up.to) {
    parent <- trunc (i/2)
    grab <- mat[parent,]
    new <- c(grab[grab != 0], i)
    mat[i,1:length(new)] <- new
}
mat
}
