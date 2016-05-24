myrpartco <- function (tree, parms)
{
    frame <- tree$frame
    node <- as.numeric(row.names(frame))
    depth <- mytree.depth(node)
    is.leaf <- (frame$var == "<leaf>")
    if (length(parms)) {
        uniform <- parms$uniform
        nspace <- parms$nspace
        minbranch <- parms$minbranch
    }
    else {
        uniform <- TRUE
        nspace <- -1
        minbranch <- 0.3
    }
    if (uniform) 
        y <- (1 + max(depth) - depth)/max(depth, 4L)
    else {
        y <- dev <- frame$dev
        temp <- split(seq(node), depth)
        parent <- match(node%/%2L, node)
        sibling <- match(ifelse(node%%2L, node - 1L, node + 1L), 
            node)
        for (i in temp[-1L]) {
            temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            y[i] <- y[parent[i]] - temp2
        }
        fudge <- minbranch * diff(range(y))/max(depth)
        for (i in temp[-1L]) {
            temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
            y[i] <- y[parent[i]] - ifelse(temp2 <= fudge & haskids, 
                fudge, temp2)
        }
        y <- y/(max(y))
    }
    x <- double(length(node))
    x[is.leaf] <- seq(sum(is.leaf))
    left.child <- match(node * 2L, node)
    right.child <- match(node * 2L + 1L, node)
    temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
    for (i in rev(temp)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
    if (nspace < 0) 
        return(list(x = x, y = y))
    compress <- function(x, me, depth) {
        lson <- me + 1L
        if (is.leaf[lson]) 
            left <- list(left = x[lson], right = x[lson], depth = depth + 
                1L, sons = lson)
        else {
            left <- compress(x, me + 1L, depth + 1L)
            x <- left$x
        }
        rson <- me + 1L + length(left$sons)
        if (is.leaf[rson]) 
            right <- list(left = x[rson], right = x[rson], depth = depth + 
                1L, sons = rson)
        else {
            right <- compress(x, rson, depth + 1L)
            x <- right$x
        }
        maxd <- max(left$depth, right$depth) - depth
        mind <- min(left$depth, right$depth) - depth
        slide <- min(right$left[1L:mind] - left$right[1L:mind]) - 
            1L
        if (slide > 0) {
            x[right$sons] <- x[right$sons] - slide
            x[me] <- (x[right$sons[1L]] + x[left$sons[1L]])/2
        }
        else slide <- 0
        if (left$depth > right$depth) {
            templ <- left$left
            tempr <- left$right
            tempr[1L:mind] <- pmax(tempr[1L:mind], right$right - 
                slide)
        }
        else {
            templ <- right$left - slide
            tempr <- right$right - slide
            templ[1L:mind] <- pmin(templ[1L:mind], left$left)
        }
        list(x = x, left = c(x[me] - nspace * (x[me] - x[lson]), 
            templ), right = c(x[me] - nspace * (x[me] - x[rson]), 
            tempr), depth = maxd + depth, sons = c(me, left$sons, 
            right$sons))
    }
    x <- compress(x, 1L, 1L)$x
    list(x = x, y = y)
}
