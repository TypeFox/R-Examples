marginMatrix <-
function (x, margin, dim = NULL, incols = FALSE, reorder = FALSE) 
{
    d = 2 - incols
    if (is.null(dim)) 
        dim = rep(2, log2(dim(x)[d]))
    if (prod(dim) != dim(x)[d]) 
        stop("Dimensions do not match")
    if (length(margin) == 0) {
        out = if (incols) 
            matrix(colSums(x), nrow = 1)
        else matrix(rowSums(x), ncol = 1)
        return(out)
    }
    else if (length(margin) == length(dim)) {
        if (reorder) {
            patt = seq_len(prod(dim))
            dim(patt) = dim
            patt = aperm(patt, rank(margin))
            dim(patt) = NULL
            if (incols) 
                x = x[patt, , drop = FALSE]
            else x = x[, patt, drop = FALSE]
        }
        return(x)
    }
    else if (length(margin) > length(dim)) 
        stop("Margin specified is too long")
    idx = array(seq_len(dim(x)[d]), dim)
    rest = seq_along(dim)[-margin]
    combs = combinations(dim[rest]) + 1
    if (reorder) {
        patt = seq_len(prod(dim[margin]))
        dim(patt) = dim[sort.int(margin)]
        patt = aperm(patt, rank(margin))
        dim(patt) = NULL
    }
    if (incols) {
        out = matrix(0, ncol = ncol(x), nrow = prod(dim[margin]))
        for (i in seq(from = 1, to = nrow(combs))) {
            init = c(subtable(idx, rest, combs[i, ]))
            out = out + x[init, , drop = FALSE]
        }
        if (reorder) 
            out = out[patt, , drop = FALSE]
    }
    else {
        out = matrix(0, nrow = nrow(x), ncol = prod(dim[margin]))
        for (i in seq(from = 1, to = nrow(combs))) {
            init = c(subtable(idx, rest, combs[i, ]))
            out = out + x[, init, drop = FALSE]
        }
        if (reorder) 
            out = out[, patt, drop = FALSE]
    }
    out
}
