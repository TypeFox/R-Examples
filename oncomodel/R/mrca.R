`mrca` <-
function (x, y, tree) 
{
    stopifnot(x %in% tree[2, ] & y %in% tree[2, ])
    if (x == y) 
        return(x)
    else {
        edge <- which(tree[2, ] == x)
        xpath <- tree[c(2, 1), edge]
        while (any(tree[2, ] == tree[1, edge])) {
            edge <- which(tree[2, ] == tree[1, edge])
            xpath <- c(xpath, tree[1, edge])
        }
        ynew <- tree[1, which(tree[2, ] == y)]
        while (!(ynew %in% xpath)) {
            ynew <- tree[1, which(tree[2, ] == ynew)]
        }
        return(ynew)
    }
}

